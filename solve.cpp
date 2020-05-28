#include "bds.h"
#include <iomanip>

using namespace std;
using namespace Eigen;

// https://twitter.com/Cshearer41/status/1265209247031271424
class CShearer20200526 : public Objective
{
public:
  double r_;
  double s_;
  double a_;
  double b_;
  double d_;
  double alpha_;
  double a0_;
  double a1_;
  double b0_;
  double b1_;
  
  CShearer20200526()
  {
  }

  void compute(const VectorXd& vars)
  {
    r_ = vars[0];
    alpha_ = vars[1];
    s_ = 1.0;
    d_ = 3*r_ + 2*s_;
    a_ = d_ * cos(alpha_);
    b_ = d_ * sin(alpha_);

    a0_ = s_ / sin(alpha_);
    a1_ = s_ * cos(alpha_);
    b0_ = s_ * sin(alpha_);
    b1_ = s_ / cos(alpha_);
  }
  
  double operator()(const VectorXd& vars)
  {
    compute(vars);
    
    double val = 0;
    val += equalityPenalty(a_*a_ + b_*b_, d_*d_);
    val += equalityPenalty(b_, s_*sin(alpha_) + s_/cos(alpha_));
    val += equalityPenalty(a_, (3*r_ + 2*s_)*cos(alpha_));
    val += equalityPenalty(a0_*a0_, s_*s_ + pow(2*r_ + s_, 2));
    val += equalityPenalty(s_*s_, a1_*a1_ + b0_*b0_);
    val += equalityPenalty(r_*r_ + s_*s_, b1_*b1_);
    return val;
  }

  string reportBest(const VectorXd& best, const string& prefix = "")
  {
    compute(best);

    double shaded = 2;
    double rect = a_*b_;
    double frac = shaded / rect;
    
    
    ostringstream oss;
    oss << prefix << "frac = " << std::setprecision(9) << frac << endl;
    oss << prefix << "a = " << a_ << endl;
    oss << prefix << "b = " << b_ << endl;
    oss << prefix << "s = " << s_ << endl;
    oss << prefix << "r = " << r_ << endl;
    oss << prefix << "d = " << d_ << endl;
    oss << prefix << "alpha = " << alpha_ << " or " << rad2deg(alpha_) << " deg" << endl;
    return oss.str();
  }
  
  VectorXd boundsLower() const
  {
    VectorXd lower(2);
    lower[0] = 0.01;  // r_
    lower[1] = 0.0;  // alpha_
    return lower;
  }
  
  VectorXd boundsUpper() const
  {
    VectorXd upper(2);
    upper[0] = 1.0;
    upper[1] = M_PI / 2.0;
    return upper;
  }

private:
};

int main(int argc, char** argv)
{
  CShearer20200526 obj;
  BigDumbSolver bds(obj, 100, 0.8);
  bds.solve();
  return 0;
}
