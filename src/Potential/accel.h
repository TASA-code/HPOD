#ifndef  _accel_H
#define  _accel_H

#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>

using namespace Eigen;

typedef Eigen::Matrix<double,6,1> Vector6d;

extern Vector3d AccelMod(Vector6d r_GCRF, double mu, int n_max, 
              int m_max, double R_ref, double time);


#endif