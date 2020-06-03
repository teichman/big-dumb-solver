#include <bds.h>
#include <iomanip>

using namespace Eigen;
using namespace std;


class OBJECTIVE_TEMPLATE : public Objective
{
public:
  // DECLARE_VARS
     
  double operator()(const VectorXd& vars)
  {
    // SET_VARS
    
    double val = 0;
    // SET_CONSTRAINTS
    return val;
  }

  string status(const string& prefix = "")
  {
    ostringstream oss;
    // PRINT_STATUS
    return oss.str();
  }
  
  VectorXd boundsLower() const
  {
    // LOWER_BOUNDS
    return lower;
  }
  
  VectorXd boundsUpper() const
  {
    // UPPER_BOUNDS
    return upper;
  }
};


int main(int argc, char** argv)
{
  OBJECTIVE_TEMPLATE obj;
  BigDumbSolver bds(obj, PARAMS);
  bds.solve();

  // TESTS
}

