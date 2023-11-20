#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>

using namespace Eigen;

typedef Matrix<double,6,1> Vector6d;

/**
* @class 
*
*/
class Coordinate {

public:

    static void check();

    static Vector6d P2ECI(Vector6d& Perifocal);
    
    static Vector3d ECI2LVLH(Vector3d& ECI_r, Vector3d& ECI_v);
    
    static double GMST(double currentTime);

    static Vector6d ECI2ECEF(Vector6d& ECI, double GMST_deg);
    
    static Vector3d SOFA_ECI2ECEF(Vector6d& ECI);
    
    static Vector2d ECEF2GEO(Vector6d& ECEF);

};
