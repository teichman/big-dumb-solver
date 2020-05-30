#pragma once

#include <Eigen/Eigen>
#include <iostream>
#include <sstream>
#include <limits>
#include <string>
#include <vector>
#include <cassert>
#include <math.h>

class Objective
{
public:
  Objective() {};
  
  virtual double operator()(const Eigen::VectorXd& vars) = 0;
  virtual Eigen::VectorXd boundsLower() const = 0;
  virtual Eigen::VectorXd boundsUpper() const = 0;
  virtual std::string status(const std::string& prefix = "") = 0;

  int dimension() const { return boundsLower().rows(); }
  std::string status(const std::string& prefix = "") const;
  void check() const;
  static double equalityPenalty(double v0, double v1) { return fabs(v0 - v1); };
};


/* Unfortunately resolution is a key sticking point.
 * Sometimes resolution of 10 just doesn't find the right answer.
 * But sometimes higher resolution is just too slow.
 */
class BigDumbSolver
{
public:
  bool quiet_;
  
  BigDumbSolver(Objective& obj, int resolution = 10, double scale_factor = 0.8, double tol = 1e-9) :
    quiet_(false),
    obj_(obj),
    resolution_(resolution),
    scale_factor_(scale_factor),
    tol_(tol)    
  {
    int num_grid_pts = pow(resolution_, obj_.dimension());
    //std::cout << "num_grid_pts: " << num_grid_pts << std::endl;
    grid_.resize(num_grid_pts, Eigen::VectorXd::Zero(obj_.dimension()));

    ticks_.resize(obj_.dimension());
    for (size_t i = 0; i < ticks_.size(); ++i)
      ticks_[i] = Eigen::VectorXd::Zero(resolution_);
  }

  Eigen::VectorXd solve();

private:
  Objective& obj_;
  int resolution_;  // Number of test points per dimension in a grid.
  double scale_factor_;  // Multiply the grid size by this amount each iteration.
  double tol_;  // Objective function must go below this to stop.

  std::vector<Eigen::VectorXd> grid_;  // Persistent storage for the grid.
  std::vector<Eigen::VectorXd> ticks_;  // Each Eigen::VectorXd is resolution_ elements of tick locations.

  void buildGrid(const Eigen::VectorXd& lower, const Eigen::VectorXd& upper,
                 int resolution, std::vector<Eigen::VectorXd>* grid_ptr,
                 std::vector<Eigen::VectorXd>* ticks_ptr) const;
  void evaluateOnImplicitGrid(const Eigen::VectorXd& lower, const Eigen::VectorXd& upper,
                              Eigen::VectorXd* best_pt, double* best_val);
};



inline double deg2rad(double deg)
{
  return deg / 180.0 * M_PI;
}

inline double rad2deg(double rad)
{
  return (rad / M_PI) * 180.0;
}

