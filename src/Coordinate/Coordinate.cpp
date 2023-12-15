#include <cmath>
#include <iostream>
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>

// #include "sofa/c/src/sofa.h"
// #include "novas/eph_manager.h"

#include "Coordinate.h"
#include "../Orbit/Orbit.h"

using namespace Eigen;  


namespace {

    VectorXi Extract_Date(std::string Date){
        std::istringstream iss(Date);

        int day, month, year, hour, minute, second;
        iss >> day;
        iss.ignore();
        std::string monthstr;
        iss >> monthstr;
        iss.ignore();
        iss >> year;
        iss.ignore();
        iss >> hour;
        iss.ignore();
        iss >> minute;
        iss.ignore();
        iss >> second;

        const char* monthCstr = monthstr.c_str();
        const char* monthName[] = {"Jan","Feb","Mar","Apr","May",
                            "Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
        
        for (int i = 0; i < 12; ++i){
            if (strcmp(monthCstr, monthName[i]) == 0){
                month = i+1;
                break;
            }
        }

        VectorXi intDate(6);
        intDate << day, month, year, hour, minute, second;

        return intDate;
    }
}



/**
* @ brief: Transform Perifocal Coordinate to ECI
* 
* @ Param: Vector6d& Perifocal, input perfifocal vector
*
*/
Vector6d Coordinate::P2ECI(Vector6d& Perifocal){
    
    Vector3d r = Perifocal.head<3>();
    Vector3d v = Perifocal.tail<3>();

    double RAAN = Orbit::RAAN;
    double i    = Orbit::i;
    double w    = Orbit::w;


    Matrix3d A_RAAN, A_i, A_w, P_ECI_Matrix;

    A_w <<  cos(w), sin(w), 0,
           -sin(w), cos(w), 0,
           0, 0, 1;

    A_i << 1, 0, 0,
           0,  cos(i), sin(i),
           0, -sin(i), cos(i);

    A_RAAN << cos(RAAN),  sin(RAAN), 0, 
              -sin(RAAN), cos(RAAN), 0, 
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




Vector6d Coordinate::ECI2ECEF(Vector6d &ECI, double t){
    
    Vector3d r = ECI.head<3>();
    Vector3d v = ECI.tail<3>();

    VectorXi Date(7);
    Date = Extract_Date(Orbit::Start_Date);

   
    int year, month, day, hour, minute, second;

    day = Date[0];
    month = Date[1];
    year = Date[2];
    hour = Date[3];
    minute = Date[4];
    second = Date[5]+t;
    // year = 2023;
    // month = 10;
    // day = 29;
    // hour = 7;
    // minute = 55;
    // second = 48+t;


    double Julian = UTC2Julian(year, month, day, hour, minute, second); 
    double theta = GMST(Julian);
    // double theta = 7.2921159e-5 * t;
    
    // Vector2d result = Precession_Nutation(Julian);
    // double deltaPsi = result[0];
    // double deltaEpsilon = result[1];

    // theta += deltaPsi * cos(deltaEpsilon);

    Matrix3d ECI_ECEF_Matrix;
    ECI_ECEF_Matrix << cos(theta), -sin(theta), 0,
                       sin(theta),  cos(theta), 0,
                       0, 0, 1;
    Vector6d ECEF;
    ECEF << ECI_ECEF_Matrix.transpose() * r, 
            ECI_ECEF_Matrix.transpose() * v;

    return ECEF;
}



Vector2d Coordinate::ECEF2GEO(Vector6d& ECEF){

    
    //ECEF.array() *= Orbit::DU;
    
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



double Coordinate::UTC2Julian(int year, int month, int day, int hour, int minute, int second){

    if (month <= 2){
        month += 12;
        year--;
    }

    // int A = year/100;
    // int B = 2 - A + A / 4;

    // Calculate Julian Date
    double JD = static_cast<int>(365.25 * year) + static_cast<int>(year/400) 
                - static_cast<int>(year/100);

    JD += static_cast<int>(30.59 * (month - 2)) + day + 1721088.5;
    JD += hour / 24.0 + minute / 1440.0 + second / 86400.0;

    double MJD = JD - 2400000.5;
    double TJD = MJD - 40000;

    //
    // double Julian = static_cast<int>(365.25 * (year + 4716)) 
    //                 + static_cast<int>(30.6001 * (month + 1));
    //
    // Julian += day + B;
    // Julian += (hour - 12) / 24.0 + minute / 1440.0 + second / 86400.0; 

    return TJD;
}




double Coordinate::GMST(double JulianDate){
    // const double T = (JulianDate - 2451545.0) / 36525.0;
    // const double theta = 280.46061837 + 360.98564736629 * (JulianDate - 2451545.0) + T * (T * (0.000387933 - T / 38710000.0));

    double theta = 280.46061837 + 360.985647 * (JulianDate - 2451545.0);

    double GMST = fmod(theta, 360.0);
    // return fmod(theta, 360.0) * M_PI / 180.0;
    return GMST * M_PI / 180;
}



Vector2d Coordinate::Precession_Nutation(double JD){

    double T = (JD - 2451545.0) / 36525.0;
    double deltaPsi_ = 5029.0966 * T + 1.11113 * T * T - 0.0000067 * T * T * T;
    
    double omega = 125.04452 - 1934.136261 * T + 0.0020708 * T * T 
                    + T * T * T / 450000.0;
    omega *= M_PI/180;
    double LO = 280.46646 + 36000.76983 * T + 0.0003032 * T * T;
    LO *= M_PI/180;
    double LS = 218.3165 + 481267.8813 * T;
    LS *= M_PI/180;

    double dPsi = -17.20 * sin(omega) - 1.32 * sin(2*LO) 
                    - 0.23 * sin(2*LS) + 0.21 * sin(2*omega);
    double dEpsilon = 9.20 * cos(omega) + 0.57 * cos(2*LO) 
                    + 0.10 * cos(2*LS) - 0.09 * cos(2*omega);

    double deltaEpsilon_ = 0.0;
    deltaPsi_ += dPsi;
    deltaPsi_ /= 3600.0;
    deltaEpsilon_ += dEpsilon / 3600.0;


        
    Vector2d result;
    result << deltaPsi_ * (M_PI/180) , deltaEpsilon_ * (M_PI/180) ;

    return result; 
}
