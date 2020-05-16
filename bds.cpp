#include <iostream>
#include <sstream>
#include <limits>
#include <string>
#include <vector>
#include <Eigen/Eigen>
#include <cassert>

using namespace std;
using namespace Eigen;

class Objective
{
public:
  Objective() {};
  
  virtual double operator()(const VectorXd& vars) const = 0;
  virtual VectorXd boundsLower() const = 0;
  virtual VectorXd boundsUpper() const = 0;

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

  VectorXd boundsLower() const
  {
    return VectorXd::Ones(4) * -10;
  }
  
  VectorXd boundsUpper() const
  {
    return VectorXd::Ones(4) * 10;
  }
};


class BigDumbSolver
{
public:
  BigDumbSolver(const Objective& obj, int resolution = 10, double scale_factor = 0.5) :
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
      if (val < best_val) {
        best_val = val;
        best_pt = grid_[i];
      }
    }

    cout << "Best point so far: " << best_pt.transpose() << endl;
    cout << "Best value so far: " << best_val << endl;
    
    // Set the new upper & lower limits.
    VectorXd range = upper - lower;
    range *= scale_factor_;
    lower = best_pt - (range / 2.0);
    upper = best_pt + (range / 2.0);

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

  // for (size_t i = 0; i < grid.size(); ++i) {
  //   cout << "Grid element " << i << " : " << grid[i].transpose() << endl;
  // }
  // cin.get();
}


int main(int argc, char** argv)
{
  ExampleObjective obj;
  cout << obj.status() << endl;
  BigDumbSolver bds(obj);
  bds.solve();
  return 0;
}

   
  
