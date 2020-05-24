#include <iostream>
#include <sstream>
#include <limits>
#include <string>
#include <vector>
#include <Eigen/Eigen>
#include <cassert>
#include <math.h>

using namespace std;
using namespace Eigen;



class Objective
{
public:
  Objective() {};
  
  virtual double operator()(const VectorXd& vars) = 0;
  virtual VectorXd boundsLower() const = 0;
  virtual VectorXd boundsUpper() const = 0;
  virtual string reportBest(const VectorXd& best, const string& prefix = "") = 0;

  int dimension() const { return boundsLower().rows(); }
  string status(const string& prefix = "") const;
  void check() const;
  static double equalityPenalty(double v0, double v1) { return fabs(v0 - v1); };
};

void Objective::check() const
{
  assert(boundsLower().rows() == boundsUpper().rows());
  assert(boundsLower().cols() == 1);
  assert(boundsUpper().cols() == 1);
  for (int i = 0; i < dimension(); ++i)
    assert(boundsUpper()(i) > boundsLower()(i));
}

string Objective::status(const string& prefix) const
{
  ostringstream oss;
  oss << prefix << "Dimension: " << dimension() << endl;
  oss << prefix << "Lower bounds: " << boundsLower().transpose() << endl;
  oss << prefix << "Upper bounds: " << boundsUpper().transpose() << endl;
  check();
  return oss.str();
}

/* Unfortunately resolution is a key sticking point.
 * Sometimes resolution of 10 just doesn't find the right answer.
 * But sometimes higher resolution is just too slow.
 */
class BigDumbSolver
{
public:
  
  BigDumbSolver(Objective& obj, int resolution = 10, double scale_factor = 0.8, double tol = 1e-9) :
    obj_(obj),
    resolution_(resolution),
    scale_factor_(scale_factor),
    tol_(tol)    
  {
    int num_grid_pts = pow(resolution_, obj_.dimension());
    cout << "num_grid_pts: " << num_grid_pts << endl;
    grid_.resize(num_grid_pts, VectorXd::Zero(obj_.dimension()));

    ticks_.resize(obj_.dimension());
    for (size_t i = 0; i < ticks_.size(); ++i)
      ticks_[i] = VectorXd::Zero(resolution_);
  }

  VectorXd solve();

private:
  Objective& obj_;
  int resolution_;  // Number of test points per dimension in a grid.
  double scale_factor_;  // Multiply the grid size by this amount each iteration.
  double tol_;  // Objective function must go below this to stop.

  vector<VectorXd> grid_;  // Persistent storage for the grid.
  vector<VectorXd> ticks_;  // Each VectorXd is resolution_ elements of tick locations.

  void buildGrid(const VectorXd& lower, const VectorXd& upper,
                 int resolution, vector<VectorXd>* grid_ptr, vector<VectorXd>* ticks_ptr) const;
  void evaluateOnImplicitGrid(const VectorXd& lower, const VectorXd& upper, VectorXd* best_pt, double* best_val);
};

VectorXd BigDumbSolver::solve()
{
  // Start at the default bounds given by the objective.
  // These variables will be changed with each iteration.
  obj_.check();
  VectorXd lower = obj_.boundsLower();
  VectorXd upper = obj_.boundsUpper();
  assert(lower.rows() == upper.rows());
  assert(lower.cols() == 1 && upper.cols() == 1);
  
  int num_evals = 0;
  double best_val = numeric_limits<double>::max();
  VectorXd best_pt = VectorXd::Ones(obj_.dimension()) * numeric_limits<double>::max();
  int iter = 0;
  while (true) {
    cout << "============================================================" << endl;
    cout << "= Iteration " << iter << endl;
    cout << "============================================================" << endl;

    evaluateOnImplicitGrid(lower, upper, &best_pt, &best_val);
    num_evals += pow(resolution_, lower.rows());
    
    // buildGrid(lower, upper, resolution_, &grid_, &ticks_);
    
    // for (size_t i = 0; i < grid_.size(); ++i) {
    //   double val = obj_(grid_[i]);
    //   ++num_evals;
    //   //cout << grid_[i].transpose() << " :: " << val << endl;
    //   if (val < best_val) {
    //     best_val = val;
    //     best_pt = grid_[i];
    //   }
    // }


    cout << "Best values so far: " << endl;
    cout << obj_.reportBest(best_pt, "  ");
    cout << "Best value so far: " << best_val << endl;
    cout << "Num evaluations so far: " << num_evals / 1e6 << "M" << endl;
    
    if (best_val < tol_) {
      cout << "Optimization complete." << endl;
      obj_(best_pt);
      break;
    }
    
    // Set the new upper & lower limits.
    VectorXd range = upper - lower;
    range *= scale_factor_;
    lower = best_pt - (range / 2.0);
    upper = best_pt + (range / 2.0);

    // If we've gone out of bounds, shift back.
    for (int i = 0; i < lower.rows(); ++i) {
      if (lower[i] < obj_.boundsLower()[i]) {
        double shift = obj_.boundsLower()[i] - lower[i];
        lower[i] += shift;
        upper[i] += shift;
        assert(upper[i] < obj_.boundsUpper()[i]);
      }
      else if (upper[i] > obj_.boundsUpper()[i]) {
        double shift = upper[i] - obj_.boundsUpper()[i];
        lower[i] -= shift;
        upper[i] -= shift;
        assert(lower[i] > obj_.boundsLower()[i]);
      }
    }
    
    iter++;
  }

  return best_pt;
}

void BigDumbSolver::buildGrid(const VectorXd& lower, const VectorXd& upper, int resolution,
                              vector<VectorXd>* grid_ptr, vector<VectorXd>* ticks_ptr) const
{
  vector<VectorXd>& ticks = *ticks_ptr;
  vector<VectorXd>& grid = *grid_ptr;
  
  int num_grid_pts = pow(resolution_, lower.rows());
  cout << "Building grid with resolution " << resolution
       << " and total num points " << num_grid_pts << endl;
  cout << "Limits:" << endl;
  cout << "  " << lower.transpose() << endl;
  cout << "  " << upper.transpose() << endl;
  
  VectorXd range = upper - lower;
  VectorXd tick_sizes = range / resolution;
  
  // Build ticks.
  assert(ticks.size() == obj_.dimension());
  for (size_t i = 0; i < ticks.size(); ++i) {
    for (int j = 0; j < resolution; ++j) {
      ticks[i][j] = lower[i] + j * tick_sizes[i];
    }
    //cout << "Dim " << i << " :: Ticks :: " << ticks[i].transpose() << endl;
  }

  // Build the grid.
  VectorXi indices = VectorXi::Zero(obj_.dimension());
  VectorXi increment_every = VectorXi::Zero(obj_.dimension());
  for (int i = 0; i < obj_.dimension(); ++i) 
    increment_every[i] = pow(resolution, i);

  // cout << "indices: " << indices.transpose() << endl;
  // cout << "increment_every: " << increment_every.transpose() << endl;

  for (size_t i = 0; i < grid.size(); ++i) {
    for (int j = 0; j < obj_.dimension(); ++j)
      if (i > 0 && i % increment_every[j] == 0)
        indices[j] = (indices[j] + 1) % resolution;
    
    //cout << "indices: " << indices.transpose() << endl;
    for (int j = 0; j < obj_.dimension(); ++j) {
      assert(grid[i].rows() == obj_.dimension());
      grid[i][j] = ticks[j][indices[j]];
    }
  }

}

void BigDumbSolver::evaluateOnImplicitGrid(const VectorXd& lower, const VectorXd& upper, VectorXd* best_pt, double* best_val)
{
  int num_grid_pts = pow(resolution_, lower.rows());

  cout << "Evaluating on grid with resolution " << resolution_
       << " and total num points " << num_grid_pts << endl;
  cout << "Limits:" << endl;
  cout << "  " << lower.transpose() << endl;
  cout << "  " << upper.transpose() << endl;
  
  VectorXd range = upper - lower;
  VectorXd tick_sizes = range / resolution_;
  
  // Build ticks.
  assert(ticks_.size() == obj_.dimension());
  for (size_t i = 0; i < ticks_.size(); ++i) {
    for (int j = 0; j < resolution_; ++j) {
      ticks_[i][j] = lower[i] + j * tick_sizes[i];
    }
    //cout << "Dim " << i << " :: Ticks :: " << ticks[i].transpose() << endl;
  }

  // Iterate through the grid.
  VectorXi indices = VectorXi::Zero(obj_.dimension());
  VectorXi increment_every = VectorXi::Zero(obj_.dimension());
  VectorXd pt = VectorXd::Zero(obj_.dimension());
  
  for (int i = 0; i < obj_.dimension(); ++i) 
    increment_every[i] = pow(resolution_, i);

  for (size_t i = 0; i < num_grid_pts; ++i) {
    for (int j = 0; j < obj_.dimension(); ++j)
      if (i > 0 && i % increment_every[j] == 0)
        indices[j] = (indices[j] + 1) % resolution_;
    
    for (int j = 0; j < obj_.dimension(); ++j) {
      pt[j] = ticks_[j][indices[j]];
    }
    
    double val = obj_(pt);
    if (val < *best_val) {
      *best_val = val;
      *best_pt = pt;
    }
    
  }
}


class ExampleObjective : public Objective
{
  double operator()(const VectorXd& vars)
  {
    double val = 0;
    for (int i = 0; i < vars.rows(); ++i) {
      val += (vars(i) - 0.0013) * (vars(i) - 0.0013);
    }
    return val;
  }

  string reportBest(const VectorXd& best, const string& prefix = "") 
  {
    ostringstream oss;
    oss << prefix << "Best: " << best.transpose() << endl;
    return oss.str();
  }

  
  VectorXd boundsLower() const
  {
    return VectorXd::Ones(4) * -10;
  }
  
  VectorXd boundsUpper() const
  {
    return VectorXd::Ones(4) * 10;
  }
};


// https://twitter.com/Cshearer41/status/1256489376881795073
class CShearer20200502 : public Objective
{
public:

  CShearer20200502()
  {
    r_ = sqrt(12.0 / M_PI);
    p1x_ = r_ * (sqrt(2.0) + 2.0) / 2.0;
    p1y_ = r_ * (sqrt(2.0) + 2.0) / 2.0;
  }
  
  double operator()(const VectorXd& vars) 
  {
    const double& a = vars[0];
    const double& b = vars[1];
    const double& r0 = vars[2];
    const double& alpha = vars[3];

    double val = 0;
    
    val += fabs(p1x_ - a + b * cos(alpha));
    val += fabs(p1y_ - b * sin(alpha));
    val += fabs(-sqrt(2.0) * a + 2.0 * (sqrt(2.0) * r_ + r0 + r_));
    val += fabs(-a + r0 + 2 * r0 * cos(alpha) + (sqrt(2.0) * r_ + r_ + r0) / sqrt(2.0));  // horiz
    val += fabs(-a + 2 * r0 + 4 * r0 * cos(alpha));  // horiz2.  redundant
    val += fabs(-a + p1y_ + p1x_ * tan(alpha) + 2.0 * r0 / cos(alpha));  // vert
    
    return val;
  }

  string reportBest(const VectorXd& best, const string& prefix = "") 
  {
    ostringstream oss;
    oss << prefix << "a = " << best[0] << endl;
    oss << prefix << "b = " << best[1] << endl;
    oss << prefix << "r0 = " << best[2] << endl;
    oss << prefix << "alpha = " << best[3] << endl;
    return oss.str();
  }
  
  VectorXd boundsLower() const
  {
    VectorXd lower(4);
    lower[0] = 0.0;
    lower[1] = 0.0;
    lower[2] = 0.0;
    lower[3] = 0.0;
    return lower;
  }
  
  VectorXd boundsUpper() const
  {
    VectorXd upper(4);
    upper[0] = 10.0;
    upper[1] = 10.0;
    upper[2] = 10.0;
    upper[3] = 2.0 * M_PI;
    return upper;
  }

private:
  double p1x_;
  double p1y_;
  double r_;
};

// https://twitter.com/Cshearer41/status/1258681340570013698
class CShearer20200508 : public Objective
{
public:

  CShearer20200508()
  {
  }
  
  double operator()(const VectorXd& vars) 
  {
    const double& r = vars[0];
    const double& h = vars[1];
    const double& alpha = vars[2];

    double val = 0;

    val += fabs(r*r + h*h - 36.0);
    val += fabs(h / r - 2.0 * h / 6.0);
    val += fabs(cos(alpha) - r / 6.0);
    
    return val;
  }

  string reportBest(const VectorXd& best, const string& prefix = "") 
  {
    ostringstream oss;
    oss << prefix << "r = " << best[0] << endl;
    oss << prefix << "h = " << best[1] << endl;
    oss << prefix << "alpha = " << best[2] << endl;
    return oss.str();
  }
  
  VectorXd boundsLower() const
  {
    VectorXd lower(3);
    lower[0] = 0.1;
    lower[1] = 0.1;
    lower[2] = 0.0;
    return lower;
  }
  
  VectorXd boundsUpper() const
  {
    VectorXd upper(3);
    upper[0] = 10.0;
    upper[1] = 10.0;
    upper[2] = 2.0 * M_PI;
    return upper;
  }

private:
};



// https://twitter.com/Cshearer41/status/1259087361645969409
class CShearer20200509 : public Objective
{
public:

  CShearer20200509()
  {
  }
  
  double operator()(const VectorXd& vars) 
  {
    const double& alpha = vars[0];
    
    double val = 0;
    val += fabs(-2.0 * sin(alpha) + 3.0 * sin(alpha) - 2.0 * cos(alpha));
    return val;
  }

  string reportBest(const VectorXd& best, const string& prefix = "") 
  {
    ostringstream oss;
    oss << prefix << "alpha = " << best[0] << endl;
    return oss.str();
  }
  
  VectorXd boundsLower() const
  {
    VectorXd lower(1);
    lower[0] = 0.0;
    return lower;
  }
  
  VectorXd boundsUpper() const
  {
    VectorXd upper(1);
    upper[0] = 2.0 * M_PI;
    return upper;
  }

private:
};


double deg2rad(double deg)
{
  return deg / 180.0 * M_PI;
}

double rad2deg(double rad)
{
  return (rad / M_PI) * 180.0;
}

// https://twitter.com/Cshearer41/status/1260847819700809728
class CShearer20200514 : public Objective
{
public:
  double cos30_;
  double sin30_;
  
  CShearer20200514()
  {
    cos30_ = cos(deg2rad(30));
    sin30_ = sin(deg2rad(30));
  }
  
  double operator()(const VectorXd& vars) 
  {
    const double& rs = vars[0];
    const double& rb = vars[1];

    double val = 0;
    val += equalityPenalty(rs + rb, cos30_);
    val += equalityPenalty(rb / cos30_, rs / (2.0 * cos30_) + sin30_);
    return val;
  }

  string reportBest(const VectorXd& best, const string& prefix = "") 
  {
    const double& rs = best[0];
    const double& rb = best[1];
    
    double a = 1.0;
    double c = rs / cos30_;
    double e = rb / cos30_;
    double f = a - c;
    double parallelogram_base = f;
    double parallelogram_height = a * cos30_;
    double shaded = 4.0 * parallelogram_base * parallelogram_height;
    double boxheight = 2.0 * parallelogram_height;
    double boxwidth = f + a + a * sin30_;
    double total = boxwidth * boxheight;
    double frac = shaded / total;
      
    ostringstream oss;
    oss << prefix << "rs = " << rs << endl;
    oss << prefix << "rb = " << rb << endl;
    oss << prefix << "a = " << a << endl;
    oss << prefix << "c = " << c << endl;
    
    oss << prefix << "frac = " << frac << endl;
    
    return oss.str();
  }
  
  VectorXd boundsLower() const
  {
    VectorXd lower(2);
    lower[0] = 0.0;
    lower[1] = 0.0;
    return lower;
  }
  
  VectorXd boundsUpper() const
  {
    VectorXd upper(2);
    upper[0] = 1.0;
    upper[1] = 1.0;
    return upper;
  }

private:
};

class CShearer20200523 : public Objective
{
public:
  double a_;
  double b_;
  double alpha_;
  double beta_;
  double c_;
  double d_;
  double area_;
  
  CShearer20200523()
  {
  }

  void compute(const VectorXd& vars)
  {
    a_ = vars[0];
    b_ = vars[1];
    alpha_ = vars[2];
    beta_ = vars[3];

    c_ = sqrt(b_*b_ - a_*a_);
    d_ = c_ / sin(alpha_);
  }
  
  double operator()(const VectorXd& vars)
  {
    compute(vars);
    
    double val = 0;
    val += equalityPenalty(a_, b_ * cos(beta_));
    val += equalityPenalty(c_, b_ * sin(beta_));
    val += equalityPenalty(tan(alpha_), c_ / (2 * a_ + c_));

    double t0 = a_ + b_ * cos(beta_);
    double t1 = b_ * sin(beta_);
    val += equalityPenalty(36.0, t0 * t0 + t1 * t1);

    val += equalityPenalty(36.0, 4*a_*a_ + c_*c_);
    
    // val += equalityPenalty(d_ * cos(alpha_), 2 * a_ + c_);
    // val += equalityPenalty(d_*d_, (2.0 * a_ + c_) * (2.0 * a_ + c_) + c_*c_);
    
    // c should be smaller than a.
    // But this probably indicates I haven't pinned down the problem yet.
    // I really need to detect multiple solutions...
    //val += max<double>(0, c_ - a_);  
    
    return val;
  }

  string reportBest(const VectorXd& best, const string& prefix = "")
  {
    compute(best);

    double gamma = atan2(b_ * sin(beta_), a_ + b_ * cos(beta_));

    // area is two squares plus two of /Delta_b minus the big long triangle on the bottom.
    area_ = (b_*b_ + a_*a_) + (c_ * a_) - (c_ * (2 * a_ + c_) / 2.0);
    
    ostringstream oss;
    oss << prefix << "a = " << a_ << endl;
    oss << prefix << "b = " << b_ << endl;
    oss << prefix << "c = " << c_ << endl;
    oss << prefix << "d = " << d_ << endl;
    oss << prefix << "alpha = " << alpha_ << " or " << rad2deg(alpha_) << " deg" << endl;
    oss << prefix << "beta = " << beta_ << " or " << rad2deg(beta_) << " deg" << endl;
    oss << prefix << "gamma = " << gamma << " or " << rad2deg(gamma) << " deg" << endl;
    oss << prefix << "area = " << area_ << endl;
    
    return oss.str();
  }
  
  VectorXd boundsLower() const
  {
    VectorXd lower(4);
    lower[0] = 0.1;  // This is what prevents it from saying the two squares are identical and tangent.
    lower[1] = 0.1;
    lower[2] = 0.0;
    lower[3] = 0.0;
    return lower;
  }
  
  VectorXd boundsUpper() const
  {
    VectorXd upper(4);
    upper[0] = 6.0;
    upper[1] = 6.0;
    upper[2] = M_PI / 2.0;
    upper[3] = M_PI / 2.0;
    return upper;
  }

private:
};

class CShearer20200430 : public Objective
{
public:
  double a_;
  double rb_;
  double rs_;
  double alpha_;
  double frac_;
  
  CShearer20200430()
  {
  }

  void compute(const VectorXd& vars)
  {
    rb_ = 1.0;
    
    a_ = vars[0];
    rs_ = vars[1];
    alpha_ = vars[2];
  }
  
  double operator()(const VectorXd& vars)
  {
    compute(vars);
    
    double val = 0;
    val += equalityPenalty(rs_*rs_ + a_*a_, 4.0 * rs_*rs_);
    val += equalityPenalty(rs_*rs_ + a_*a_ / 4.0, 1.0);
    val += equalityPenalty(sin(alpha_), rs_);
    val += equalityPenalty(cos(alpha_), a_ / 2.0);
    return val;
  }

  string reportBest(const VectorXd& best, const string& prefix = "")
  {
    compute(best);

    frac_ = rs_*rs_;
    
    ostringstream oss;
    oss << prefix << "a = " << a_ << endl;
    oss << prefix << "rs = " << rs_ << endl;
    oss << prefix << "rb = " << rb_ << endl;
    oss << prefix << "alpha = " << alpha_ << " or " << rad2deg(alpha_) << " deg" << endl;
    oss << prefix << "frac = " << frac_ << endl;
    return oss.str();
  }
  
  VectorXd boundsLower() const
  {
    VectorXd lower(3);
    lower[0] = 0.0;
    lower[1] = 0.0;
    lower[2] = 0.0;
    return lower;
  }
  
  VectorXd boundsUpper() const
  {
    VectorXd upper(3);
    upper[0] = 2.0;
    upper[1] = 1.0;
    upper[2] = M_PI / 2.0;
    return upper;
  }

private:
};

// TODO: Proper test infrastructure.
void test()
{
  double tol = 1e-6;
  
  { 
    CShearer20200523 obj;
    BigDumbSolver bds(obj, 10, 0.5);
    bds.solve();
    assert(fabs(obj.area_ - 18.0) < tol);
  }

  {
    CShearer20200430 obj;
    BigDumbSolver bds(obj, 30, 0.75);
    bds.solve();
    assert(fabs(obj.frac_ - 4.0/7.0) < tol);
  }
  
  cout << "Tests complete." << endl;
  exit(0);
}


int main(int argc, char** argv)
{
  test();

  // ExampleObjective obj;
  //CShearer20200502 obj;
  //CShearer20200508 obj;
  //CShearer20200509 obj;
  //CShearer20200514 obj;
  //CShearer20200523 obj;
  CShearer20200430 obj;
  cout << obj.status() << endl;
  
  BigDumbSolver bds(obj, 30, 0.75);
  bds.solve();
  return 0;
}

   
  
