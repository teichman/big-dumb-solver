#include "bds.h"

using namespace std;
using namespace Eigen;

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
    if (!quiet_) {
      cout << "============================================================" << endl;
      cout << "= Iteration " << iter << endl;
      cout << "============================================================" << endl;
    }

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
    
    if (!quiet_) {
      cout << "Best values so far: " << endl;
      cout << obj_.reportBest(best_pt, "  ");  
      cout << "Best value so far: " << best_val << endl;
      cout << "Num evaluations so far: " << num_evals / 1e6 << "M" << endl;
    }
    
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
  if (!quiet_) {
    cout << "Building grid with resolution " << resolution
         << " and total num points " << num_grid_pts << endl;
    cout << "Limits:" << endl;
    cout << "  " << lower.transpose() << endl;
    cout << "  " << upper.transpose() << endl;
  }
  
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

  if (!quiet_) {
    cout << "Evaluating on grid with resolution " << resolution_
         << " and total num points " << num_grid_pts << endl;
    cout << "Limits:" << endl;
    cout << "  " << lower.transpose() << endl;
    cout << "  " << upper.transpose() << endl;
  }
  
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
