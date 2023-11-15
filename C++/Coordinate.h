#ifndef COORINDATE_H
#define COORDINATE_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>


/**
* @class 
*
*/
class Coordinate {
    


public:

    Eigen::Vector3d LVLH_r;
    Eigen::Vector3d LVLH_v;

    Eigen::Matrix3d ECI_LVLH;
        
    Eigen::Vector3d P2ECI();

    Eigen::Vector3d ECI2ECEF();

    Eigen::Vector3d ECI2LVLH(Eigen::Vector3d&ECI_r, Eigen::Vector3d&ECI_v);


};

#endif // COORINDATE_H
