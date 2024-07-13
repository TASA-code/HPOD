#ifndef  _rk45_H
#define  _rk45_H

#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>
using namespace Eigen;
typedef Eigen::Matrix<double,6,1> Vector6d;


extern void RungeKutta45(const double& T, const double& dt, const int& outputFrequency, Vector6d& x);



#endif