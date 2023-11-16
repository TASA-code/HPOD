#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>

#include "Orbit.h"
using namespace Eigen;

typedef Matrix<double,6,1> Vector6d;

/**
* @class 
*
*/
class Coordinate : public Orbit {

public:

    Vector3d LVLH_r;
    Vector3d LVLH_v;

    Matrix3d ECI_LVLH;
        
    Vector3d P2ECI();
    
    Vector3d ECI2LVLH(Vector3d& ECI_r, Vector3d& ECI_v);
    
    Vector3d ECI2ECEF(Vector6d& ECI, double GMST_deg);

    Vector3d ECEF2GEO(Vector3d& ECEF);

};
