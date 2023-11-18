#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>

#include "Coordinate.h"
#include "Orbit.h"


using namespace Eigen;



void Coordinate::check(){
    std::cout << "checkpoint" << std::endl;
}


Vector6d Coordinate::P2ECI(Vector6d& Perifocal){
    
    Vector3d r = Perifocal.head<3>();
    Vector3d v = Perifocal.tail<3>();

    double RAAN = Orbit::RAAN;
    double i    = Orbit::i;
    double w    = Orbit::w;


    Matrix3d A_RAAN, A_i, A_w, P_ECI_Matrix;

    A_RAAN << cos(RAAN),  sin(RAAN), 0, 
              -sin(RAAN), cos(RAAN), 0, 
              0, 0, 1;

    A_i << 1, 0, 0,
           0,  cos(i), sin(i),
           0, -sin(i), cos(i);

    A_w <<  cos(w), sin(w), 0,
           -sin(w), cos(w), 0,
           0, 0, 1;

    P_ECI_Matrix = (A_w * A_i * A_RAAN).transpose();

    Vector3d ECI_r, ECI_v;
    ECI_r = P_ECI_Matrix * r;
    ECI_v = P_ECI_Matrix * v;
        
    Vector6d ECI;
    ECI << ECI_r, ECI_v;

    return ECI;
 }




Vector3d Coordinate::ECI2LVLH(Vector3d&ECI_r, Vector3d&ECI_v){
    
    Vector3d LV = ECI_r.normalized();

    Vector3d h = LV.cross(ECI_v);
    
    Vector3d orbit_normal = h.normalized();
    
    Vector3d LH = orbit_normal.cross(LV);

    Matrix3d ECI_LVLH_Matrix; 
    ECI_LVLH_Matrix << LV, LH, orbit_normal;

    Vector3d LVLH_r, LVLH_v;
    LVLH_r = ECI_LVLH_Matrix * ECI_r;
    LVLH_v = ECI_LVLH_Matrix * ECI_v;
    
    return LVLH_r;
}


//
// std::vector<double> Orbit::ECI2LVLH(std::vector<double>&ECI){
//     
//     std::vector<double> ECI_r = {0.0,0.0,0.0};
//     std::vector<double> ECI_v = {0.0,0.0,0.0};
//     for(int i = 0; i < 3; i++){
//         ECI_r[i] = ECI[i];
//         ECI_v[i] = ECI[i+3];
//     }
//
//     double r_norm = sqrt(ECI_r[0] * ECI_r[0] + ECI_r[1] * ECI_r[1] + ECI_r[2] * ECI_r[2]);
//     std::vector<double> LV = {0.0, 0.0, 0.0};
//     
//     for (int i = 0; i < 3; ++i){
//         LV[i] = ECI_r[i]/r_norm;
//     };
//
//     
//     std::vector<double> h(3,0.0);
//     h = {(ECI_r[1]*ECI_v[2]-ECI_r[2]*ECI_v[1]), (ECI_r[0]*ECI_v[2]-ECI_r[2]*ECI_v[0]), (ECI_r[0]*ECI_v[1]-ECI_r[1]*ECI_v[0])};
//
//     
//     double h_norm = sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]);
//     std::vector<double> orbit_normal = {0.0, 0.0, 0.0};
//     for (int i = 0; i<3; ++i){
//         orbit_normal[i] = h[i]/h_norm;
//     };
//
//
//     std::vector<double> LH(3,0.0); 
//     LH = {(orbit_normal[1]*LV[2]-orbit_normal[2]*LV[1]), (orbit_normal[0]*LV[2]-orbit_normal[2]*LV[0]), (orbit_normal[0]*LV[1]-orbit_normal[1]*LV[0])};
//
//     std::vector<std::vector<double>> ECILVLH = {
//                 {LV[0], LH[0], orbit_normal[0]},
//                 {LV[1], LH[1], orbit_normal[1]},
//                 {LV[2], LH[2], orbit_normal[2]},
//                                                 };
//     std::vector<double> LVLH_x = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//
//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++){
//             LVLH_x[i] += ECILVLH[i][j] * ECI_r[j];
//             LVLH_x[i+3] += ECILVLH[i][j] * ECI_v[j];
//         }
//     }
//
//     return LVLH_x;
// }
//
// std::vector<double> quaternion(std::vector<double> LVLH){
//
//     std::vector<double> r = {LVLH[0], LVLH[1], LVLH[2]};
//     double norm = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
//
//     std::vector<double> LVLH_normalised = {0.0, 0.0, 0.0};
//     for (int i = 0; i < 3; i++){
//         LVLH_normalised[i] = r[i]/norm;
//     };
//     
//     // pitch
//     double phi = -asin(LVLH_normalised[0]);
//     
//     // roll
//     double theta = atan2(LVLH_normalised[1], LVLH_normalised[2]);
//     
//     // yaw
//     double psi = atan2(LVLH_normalised[1], LVLH_normalised[0]);
//     
//
//     // quaternion
//     double qw = cos(phi/2.0);
//     double qx = sin(phi/2.0) * cos(theta/2.0);
//     double qy = sin(phi/2.0) * sin(theta/2.0) * cos(psi/2.0);
//     double qz = sin(phi/2.0) * sin(theta/2.0) * sin(psi/2.0);
//     
//     std::vector<double> quaternion = {qw,qx,qy,qz};
//
//     return quaternion;
// }



double Coordinate::GMST(double currentTime){
    
    double year, month, day, hour, minute, second;

    year = 2023.0;
    month = 10.0;
    day = 29.0;
    hour = 7;
    minute = 55;
    second = 48;

    if (month <= 2) {
       year--;
       month += 12;
     }

    // Calculate Julian Date
    int A = year / 100;
    int B = 2 - A + year;
    double julianDate = floor(365.25 * B) + floor(30.6001 * (month + 1)) + day + (hour + minute / 60.0 + (second + currentTime) / 3600.0) + 1720996.5;
    
    double T = (julianDate - 2451545.0) / 36525.0;

    // Mean sidereal time at the Greenwich meridian at J2000.0 (in seconds)
    double GMST = 280.46061837 + 360.98564736629 * (julianDate - 2451545.0) + 0.000387933 * T * T - T * T * T / 38710000.0;

    // Ensure the result is within the range [0, 360)
    GMST = fmod(GMST, 360.0);

    if (GMST < 0.0) {
        GMST += 360.0;
    }       

    return GMST * (M_PI/180); 
}



Vector6d Coordinate::ECI2ECEF(Vector6d &ECI, double GMST){
    
    Vector3d r = ECI.head<3>();
    Vector3d v = ECI.tail<3>();

    Matrix3d ECI_ECEF_Matrix;
    ECI_ECEF_Matrix << cos(GMST), -sin(GMST), 0,
                       sin(GMST),  cos(GMST), 0,
                       0, 0, 1;
    Vector6d ECEF;
    ECEF << ECI_ECEF_Matrix * r, ECI_ECEF_Matrix * v;

    return ECEF;
}



Vector2d Coordinate::ECEF2GEO(Vector6d& ECEF){

    
    ECEF.array() *= Orbit::DU;
    
    const double a = 6378.137;
    const double b = 6356.7534;
    const double e2 = 1.0 - (b*b) / (a*a);

    // Calculate geodetic latitude and longitude
    double longitude = atan2(ECEF[1], ECEF[0]) * (180.0 / M_PI);
    double p = sqrt(ECEF[0] * ECEF[0] + ECEF[1] * ECEF[1]);
    double latitude = atan2(ECEF[2], p);

    // Iterate to calculate altitude
    double latitude_prev = 0.0;
    double altitude = 0.0;

    const double epsilon = 1e-9; // Tolerance for convergence
    while (std::abs(latitude - latitude_prev) > epsilon) {
        latitude_prev = latitude;
        double N = a / sqrt(1.0 - e2 * sin(latitude) * sin(latitude));
        altitude = p / cos(latitude) - N;
        latitude = atan2(ECEF[2], p / (1.0 - e2 * N / (N + altitude)));
    }

    // Convert latitude and longitude to degrees
    latitude *= (180.0 / M_PI);
     
    Vector2d GEO;
    GEO << longitude, latitude;

    return GEO;
}





