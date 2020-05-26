#include "bds.h"

using namespace std;
using namespace Eigen;

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


int main(int argc, char** argv)
{
  CShearer20200525 obj;
  BigDumbSolver bds(obj, 30, 0.75);
  bds.solve();
  return 0;
}
