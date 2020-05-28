#include "bds.h"
#include <iomanip>

using namespace std;
using namespace Eigen;


int main(int argc, char** argv)
{
  CShearer20200526 obj;
  BigDumbSolver bds(obj, 100, 0.8);
  bds.solve();
  return 0;
}
