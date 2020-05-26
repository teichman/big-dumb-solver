#include "bds.h"

using namespace std;
using namespace Eigen;

class ExampleObjective : public Objective
{
  double operator()(const VectorXd& vars)
  {
    double val = 0;
    for (int i = 0; i < vars.rows(); ++i) {
      val += (vars(i) - 0.4213) * (vars(i) - 0.4213);
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
    return VectorXd::Ones(4) * -1;
  }
  
  VectorXd boundsUpper() const
  {
    return VectorXd::Ones(4) * 1;
  }
};


int main(int argc, char** argv)
{
  //test();
  ExampleObjective obj;
  //CShearer20200502 obj;
  //CShearer20200508 obj;
  //CShearer20200509 obj;
  //CShearer20200514 obj;
  //CShearer20200523 obj;
  //CShearer20200430 obj;
  cout << obj.status() << endl;
  
  BigDumbSolver bds(obj, 30, 0.75);
  bds.solve();
  return 0;
}
