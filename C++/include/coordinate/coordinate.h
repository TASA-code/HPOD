#ifndef  _coordinate_H
#define  _coordinate_H

#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>

using namespace Eigen;

typedef Matrix<double,6,1> Vector6d;

extern Vector6d OE2ECI(const double* OE);

extern Vector6d P2ECI(Vector6d& Perifocal);

extern Vector6d ECI2ECEF(const Vector6d& ECI, double t);

extern Vector3d ECEF2ECI(Vector3d& a, double t);

extern Vector2d ECEF2GEO(Vector6d& ECEF);


#endif