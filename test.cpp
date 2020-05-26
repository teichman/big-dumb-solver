#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "bds.h"

using namespace Eigen;
using namespace std;

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


TEST_CASE("doctest::Approx")
{
  // This passes.  1e-5 fails.
  CHECK(1e-4 != doctest::Approx(0.0));
}


TEST_CASE("CShearer20200523")
{
  CShearer20200523 obj;
  BigDumbSolver bds(obj, 10, 0.5);
  bds.quiet_ = true;
  VectorXd pt = bds.solve();
  cout << "Best: " << endl;
  cout << obj.reportBest(pt, "  ");
  CHECK(obj.area_ == doctest::Approx(18.0));
}


TEST_CASE("CShearer20200430")
{
  CShearer20200430 obj;
  BigDumbSolver bds(obj, 30, 0.75);
  bds.quiet_ = true;
  VectorXd pt = bds.solve();
  cout << "Best: " << endl;
  cout << obj.reportBest(pt, "  ");
  CHECK(obj.frac_ == doctest::Approx(4.0/7.0));
}


// https://twitter.com/Cshearer41/status/1264912625378504709
class CShearer20200525 : public Objective
{
public:
  double rs_;
  double ru_;
  double rl_;
  double alpha_;
  double ay_;
  double ab_;
  
  CShearer20200525()
  {
  }

  void compute(const VectorXd& vars)
  {
    ru_ = vars[0];
    alpha_ = vars[1];
    rs_ = ru_ * cos(alpha_) / 2.0;
    rl_ = ru_ * sin(alpha_);
    ay_ = M_PI * rl_*rl_ / 2.0 - ru_*ru_*alpha_ + 2.0*rl_*rs_;
    ab_ = M_PI * ru_*ru_ / 2.0 - ru_*ru_*alpha_ + 2.0*rl_*rs_ - 2.0 * M_PI * rs_*rs_;
  }
  
  double operator()(const VectorXd& vars)
  {
    compute(vars);
    
    double val = 0;
    val += equalityPenalty(ay_, 16.0);
    val += equalityPenalty(ru_*ru_, 4.0 * rs_*rs_ + rl_*rl_);    
    return val;
  }

  string reportBest(const VectorXd& best, const string& prefix = "")
  {
    compute(best);
    
    ostringstream oss;
    oss << prefix << "ru = " << ru_ << endl;
    oss << prefix << "rs = " << rs_ << endl;
    oss << prefix << "rl = " << rl_ << endl;
    oss << prefix << "alpha = " << alpha_ << " or " << rad2deg(alpha_) << " deg" << endl;
    oss << prefix << "ay = " << ay_ << endl;
    oss << prefix << "ab = " << ab_ << endl;
    return oss.str();
  }
  
  VectorXd boundsLower() const
  {
    VectorXd lower(2);
    lower[0] = 0.0;  // ru_
    lower[1] = 0.0;  // alpha_
    return lower;
  }
  
  VectorXd boundsUpper() const
  {
    VectorXd upper(2);
    upper[0] = 10.0;
    upper[1] = M_PI / 2.0;
    return upper;
  }

private:
};

TEST_CASE("CShearer20200525")
{
  CShearer20200525 obj;
  BigDumbSolver bds(obj, 30, 0.75);
  bds.quiet_ = true;
  VectorXd pt = bds.solve();
  // cout << "Best: " << endl;
  // cout << obj.reportBest(pt, "  ");
  CHECK(obj.ab_ == doctest::Approx(16.0));
}
