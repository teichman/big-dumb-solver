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
  
  virtual double operator()(const VectorXd& vars) const = 0;
  virtual VectorXd boundsLower() const = 0;
  virtual VectorXd boundsUpper() const = 0;
  virtual string reportBest(const VectorXd& best, const string& prefix = "") const = 0;

  int dimension() const { return boundsLower().rows(); };
  string status(const string& prefix = "") const;
  void check() const;
  
private:
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

class BigDumbSolver
{
public:
  BigDumbSolver(const Objective& obj, int resolution = 10, double scale_factor = 0.8) :
    obj_(obj),
    resolution_(resolution),
    scale_factor_(scale_factor)
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
  const Objective& obj_;
  int resolution_;  // Number of test points per dimension in a grid.
  double scale_factor_;  // Multiply the grid size by this amount each iteration.

  vector<VectorXd> grid_;  // Persistent storage for the grid.
  vector<VectorXd> ticks_;  // Each VectorXd is resolution_ elements of tick locations.

  void buildGrid(const VectorXd& lower, const VectorXd& upper,
                 int resolution, vector<VectorXd>* grid_ptr, vector<VectorXd>* ticks_ptr) const;

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

  double best_val = numeric_limits<double>::max();
  VectorXd best_pt = VectorXd::Ones(obj_.dimension()) * numeric_limits<double>::max();
  int iter = 0;
  while (true) {
    cout << "============================================================" << endl;
    cout << "= Iteration " << iter << endl;
    cout << "============================================================" << endl;

    buildGrid(lower, upper, resolution_, &grid_, &ticks_);
    
    for (size_t i = 0; i < grid_.size(); ++i) {
      double val = obj_(grid_[i]);
      //cout << grid_[i].transpose() << " :: " << val << endl;
      if (val < best_val) {
        best_val = val;
        best_pt = grid_[i];
      }
    }

    cout << "Best value so far: " << best_val << endl;
    cout << "Best point so far: " << best_pt.transpose() << endl;
    cout << obj_.reportBest(best_pt, "  ");
    
    // Set the new upper & lower limits.
    VectorXd range = upper - lower;
    range *= scale_factor_;
    lower = best_pt - (range / 2.0);
    lower = lower.cwiseMax(obj_.boundsLower());
    upper = best_pt + (range / 2.0);
    upper = upper.cwiseMin(obj_.boundsUpper());

    iter++;
  }

  return VectorXd::Zero(obj_.dimension());
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
    cout << "Dim " << i << " :: Ticks :: " << ticks[i].transpose() << endl;
  }

  // Build the grid.
  VectorXi indices = VectorXi::Zero(obj_.dimension());
  VectorXi increment_every = VectorXi::Zero(obj_.dimension());
  for (int i = 0; i < obj_.dimension(); ++i) 
    increment_every[i] = pow(resolution, i);

  cout << "indices: " << indices.transpose() << endl;
  cout << "increment_every: " << increment_every.transpose() << endl;

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


class ExampleObjective : public Objective
{
  double operator()(const VectorXd& vars) const
  {
    double val = 0;
    for (int i = 0; i < vars.rows(); ++i) {
      val += (vars(i) - 0.0013) * (vars(i) - 0.0013);
    }
    return val;
  }

  string reportBest(const VectorXd& best, const string& prefix = "") const
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


class CShearer20200502Obj : public Objective
{
public:

  CShearer20200502Obj()
  {
    r_ = sqrt(12.0 / M_PI);
    p1x_ = r_ * (sqrt(2.0) + 2.0) / 2.0;
    p1y_ = r_ * (sqrt(2.0) + 2.0) / 2.0;
  }
  
  double operator()(const VectorXd& vars) const
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
    //val += fabs(-a + p1y_ + p1x_ * tan(alpha) + 2.0 * r0 / cos(alpha));  // vert
    
    return val;
  }

  string reportBest(const VectorXd& best, const string& prefix = "") const
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
    upper[0] = 100.0;
    upper[1] = 100.0;
    upper[2] = 100.0;
    upper[3] = 2.0 * M_PI;
    return upper;
  }

private:
  double p1x_;
  double p1y_;
  double r_;
};

// https://twitter.com/Cshearer41/status/1258681340570013698
class CShearer20200508Obj : public Objective
{
public:

  CShearer20200508Obj()
  {
  }
  
  double operator()(const VectorXd& vars) const
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

  string reportBest(const VectorXd& best, const string& prefix = "") const
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




int main(int argc, char** argv)
{
  // ExampleObjective obj;
  // cout << obj.status() << endl;

  //CShearer20200502Obj obj;
  CShearer20200508Obj obj;
  cout << obj.status() << endl;
  
  BigDumbSolver bds(obj, 100);
  bds.solve();
  return 0;
}

   
  
