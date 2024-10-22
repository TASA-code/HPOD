#include <cmath>
#include <iostream>

#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>

#include "coordinate.h"
#include "propagator.h"
#include "dm_time.h"

using namespace Eigen;  

typedef Eigen::Matrix<double,6,1> Vector6d;

// Helper function to unpack Eigen::Map into tuple
template <typename Derived>
auto unpackMap(const Eigen::DenseBase<Derived>& map) {
    return std::make_tuple(map[0], map[1], map[2], map[3], map[4], map[5]);
}




double GMST(double JulianDate){

    /* Local variables */
    double theta;           /* Rotation Angle */
    double GMST;            /* Greenwich Mean Sidereal Time */


    theta = 280.46061837 + 360.985647 * (JulianDate - 2451545.0);
    GMST = fmod(theta, 360.0);

    return GMST * M_PI / 180;
}



double UTC2Julian(int year, int month, int day, int hour, int minute, int second){

    double JD;              /* Julian Day */
    double MJD;             /* Modified Julian Day */
    double TJD;             /* Truncated Julian Day */


    if (month <= 2){
        month += 12;
        year--;
    }

    // Calculate Julian Date
    JD = static_cast<int>(365.25 * year) + static_cast<int>(year/400) 
                - static_cast<int>(year/100);

    JD += static_cast<int>(30.59 * (month - 2)) + day + 1721088.5;
    JD += hour / 24.0 + minute / 1440.0 + second / 86400.0;

    MJD = JD - 2400000.5;
    TJD = MJD - 40000;

    return TJD;
}





/**
* @ brief: Transform Perifocal Coordinate to ECI
* 
* @ Param: Vector6d& Perifocal, input perfifocal vector
*
*/
Vector6d OE2ECI(const double* OE){

    const double Earth_mu = Propagator::Earth_mu;
    Vector6d state;
    // Define Eigen::Map for direct access to OE
    Map<const VectorXd> variables(OE, 6); // Map OE to a const Eigen vector of size 6

    // Unpack variables using helper function and std::tie
    double SMA, e, i, M, w, RAAN;
    std::tie(SMA, e, i, M, w, RAAN) = unpackMap(variables);


    // Iteration to convert mean anomaly to eccentric anomaly
    double E = M;
    double f = E - e * sin(E) - M;

    while (fabs(f) > 1e-9){ 
        E = E - f / (1.0 - e*cos(E));
        f = E - e*sin(E) - M;
    }
    
    // Convert Eccentric anomaly to true anomaly
    double theta = 2.0 * atan2(sqrt(1.0 + e) * sin(E / 2.0),
                               sqrt(1.0 - e) * cos(E / 2.0));
    
    
    // initialise vector
    Eigen::Vector3d P_r, P_v, r_vector, v_vector;
    

    double h = sqrt(Earth_mu*SMA*1000*(1-pow(e, 2.0)));

    // Define unit vectors
    Vector3d i_e(1.0, 0.0, 0.0); // i_e = [1; 0; 0]
    Vector3d i_p(0.0, 1.0, 0.0); // i_p = [0; 1; 0]

    // Calculate r0_p (position in perifocal frame)
    Vector3d r0_p = (pow(h, 2.0) / (Earth_mu * (1 + e * cos(theta)))) * (cos(theta) * i_e + sin(theta) * i_p);
 
    // Calculate v0_p (velocity in perifocal frame)
    Vector3d v0_p = (Earth_mu / h) * (-sin(theta) * i_e + (e + cos(theta)) * i_p);


    Vector6d Perifocal;
    Perifocal << r0_p, v0_p;
    state = P2ECI(Perifocal);

    return state;

};



/**
* @ brief: Transform Perifocal Coordinate to ECI
* 
* @ Param: Vector6d& Perifocal, input perfifocal vector
*
*/
Vector6d P2ECI(Vector6d& Perifocal){

    /* Local variables */
    double   RAAN;                  /* Right Ascension Ascending Node */
    double   i;                     /* Inclination */
    double   w;                     /* Argument of perigee */
    Vector3d r;                     /* Perifocal position */
    Vector3d v;                     /* Perifocal velocity*/
    Vector6d ECI;                   /* Earth-Centred inertia */
    Vector3d ECI_r, ECI_v;          /* Earth-Centred position & velocity*/
    Matrix3d A_RAAN, A_i, A_w;      /* Rotational matrix*/
    Matrix3d P_ECI_Matrix;
    /* end of Local variables */


    r = Perifocal.head<3>();
    v = Perifocal.tail<3>();
    RAAN = Propagator::RAAN;
    i    = Propagator::i;
    w    = Propagator::w;


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

    
    ECI_r = P_ECI_Matrix * r;
    ECI_v = P_ECI_Matrix * v;
    ECI << ECI_r, ECI_v;

    return ECI;
}




Vector6d ECI2ECEF(const Vector6d& ECI, double t){
    
    Vector3d r = ECI.head<3>();
    Vector3d v = ECI.tail<3>();

    VectorXi Date(7);
    Date = Extract_Date(Propagator::Start_Date);

    int year, month, day, hour, minute, second;

    day = Date[0];
    month = Date[1];
    year = Date[2];
    hour = Date[3];
    minute = Date[4];
    second = Date[5]+t;

    double Julian = UTC2Julian(year, month, day, hour, minute, second); 
    double theta = GMST(Julian);

    Matrix3d ECI_ECEF_Matrix;
    ECI_ECEF_Matrix << cos(theta), -sin(theta), 0,
                       sin(theta),  cos(theta), 0,
                       0, 0, 1;
    Vector6d ECEF;
    ECEF << ECI_ECEF_Matrix.transpose() * r, 
            ECI_ECEF_Matrix.transpose() * v;

    return ECEF;
}


Vector3d ECEF2ECI(Vector3d& a, double t){

    VectorXi Date(7);
    Date = Extract_Date(Propagator::Start_Date);

    int year, month, day, hour, minute, second;

    day = Date[0];
    month = Date[1];
    year = Date[2];
    hour = Date[3];
    minute = Date[4];
    second = Date[5]+t;


    double Julian = UTC2Julian(year, month, day, hour, minute, second); 
    double theta = GMST(Julian);

    Matrix3d ECEF_ECI_Matrix;
    ECEF_ECI_Matrix << cos(theta), -sin(theta), 0,
                       sin(theta),  cos(theta), 0,
                       0, 0, 1;

    Vector3d ECI_a;
    ECI_a << ECEF_ECI_Matrix * a;

    return ECI_a;
}


Vector2d ECEF2GEO(Vector6d& ECEF){

    const double a = 6378.137e3;
    const double b = 6356.7534e3;
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



// Vector2d Precession_Nutation(double JD){

//     double T = (JD - 2451545.0) / 36525.0;
//     double deltaPsi_ = 5029.0966 * T + 1.11113 * T * T - 0.0000067 * T * T * T;
    
//     double omega = 125.04452 - 1934.136261 * T + 0.0020708 * T * T 
//                     + T * T * T / 450000.0;
//     omega *= M_PI/180;
//     double LO = 280.46646 + 36000.76983 * T + 0.0003032 * T * T;
//     LO *= M_PI/180;
//     double LS = 218.3165 + 481267.8813 * T;
//     LS *= M_PI/180;

//     double dPsi = -17.20 * sin(omega) - 1.32 * sin(2*LO) 
//                     - 0.23 * sin(2*LS) + 0.21 * sin(2*omega);
//     double dEpsilon = 9.20 * cos(omega) + 0.57 * cos(2*LO) 
//                     + 0.10 * cos(2*LS) - 0.09 * cos(2*omega);

//     double deltaEpsilon_ = 0.0;
//     deltaPsi_ += dPsi;
//     deltaPsi_ /= 3600.0;
//     deltaEpsilon_ += dEpsilon / 3600.0;


        
//     Vector2d result;
//     result << deltaPsi_ * (M_PI/180) , deltaEpsilon_ * (M_PI/180) ;

//     return result; 
// }
