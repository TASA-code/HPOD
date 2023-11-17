#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>       
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>

#include "Orbit.h"
#include "Coordinate.h"


// static function initialisation
double Orbit::SMA   = 0.0;
double Orbit::e     = 0.0; 
double Orbit::i     = 0.0; 
double Orbit::RAAN  = 0.0;
double Orbit::w     = 0.0;
double Orbit::theta = 0.0;


/**
 * @brief Takes user input orbital element and set initial coordinates
 *
 * @param arg_SMA
 * @param arg_e
 * @param arg_i
 * @param arg_RAAN
 * @param arg_w
 * @param arg_theta
 */
void Orbit::SetParameter(const double &arg_SMA, const double &arg_e,
                         const double &arg_i, const double &arg_RAAN,
                         const double &arg_w, const double &arg_theta) {

    // Printing values of all the parameters to the terminal.
    std::cout << std::endl;
    std::cout << "\t\tList of Commands Line Inputs:\n" << std::endl;
    std::cout << "\t\t- SMA       :  " << arg_SMA << std::endl;
    std::cout << "\t\t- e         :  " << arg_e << std::endl;
    std::cout << "\t\t- i         :  " << arg_i << std::endl;
    std::cout << "\t\t- RAAN      :  " << arg_RAAN << std::endl;
    std::cout << "\t\t- w         :  " << arg_w << std::endl;
    std::cout << "\t\t- theta     :  " << arg_theta << std::endl;
    std::cout << std::endl;
    

    SMA = arg_SMA;
    e = arg_e;
    i = arg_i * M_PI / 180;
    RAAN = arg_RAAN * M_PI / 180;
    w = arg_w * M_PI / 180;
    theta = arg_theta * M_PI / 180;

    // initialise vector
    r_perifocal = {0.0, 0.0, 0.0};
    v_perifocal = {0.0, 0.0, 0.0};
    std::vector<double> i_e = {1.0, 0.0, 0.0};
    std::vector<double> i_p = {0.0, 1.0, 0.0};

    // calculate angular momentum h
    h = std::sqrt(Earth_mu * SMA * (1 - e * e));

    for (int i = 0; i < 3; ++i) {
      r_perifocal[i] = ((h * h) / (Earth_mu * (1 + e * std::cos(theta)))) *
                      (std::cos(theta) * i_e[i] + std::sin(theta) * i_p[i]);
      v_perifocal[i] = (Earth_mu / h) * (-std::sin(theta) * i_e[i] +
                                        (e + std::cos(theta)) * i_p[i]);
    }
};



/**
 * @brief Perform coordinate transform from Perifocal frame to ECI
 *
 */
void Orbit::P2ECI() {
    Coordinate::check();
    
    // Define transformation matrices R_RAAN, R_i, and R_w
    std::vector<std::vector<double>> R_RAAN = {
        {cos(RAAN), sin(RAAN), 0}, {-sin(RAAN), cos(RAAN), 0}, {0, 0, 1}};

    std::vector<std::vector<double>> R_i = {
        {1, 0, 0}, {0, cos(i), sin(i)}, {0, -sin(i), cos(i)}};

    std::vector<std::vector<double>> R_w = {
        {cos(w), sin(w), 0}, {-sin(w), cos(w), 0}, {0, 0, 1}};

    // Calculate the transformation matrix A = R_RAAN * R_i * R_w
    // Initialise matrix A
    RotMatrix_P2ECI = std::vector<std::vector<double>>(3, std::vector<double>(3, 0.0));
    RotMatrix_P2ECI = matrixmultiply(matrixmultiply(R_w, R_i), R_RAAN);

    std::vector<std::vector<double>> R_transpose(3, std::vector<double>(3, 0.0));
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        R_transpose[i][j] = RotMatrix_P2ECI[j][i];
      }
    }

    r_ECI = std::vector<double>(3, 0.0);
    v_ECI = std::vector<double>(3, 0.0);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        r_ECI[i] += R_transpose[i][j] * r_perifocal[j];
        v_ECI[i] += R_transpose[i][j] * v_perifocal[j];
      }
    }

    // Print the resulting vectors r_ECI and v_ECI
    // format data output
    std::cout << "\t\t----------------------------------" << std::endl;
    std::cout << "\t\t...PRINTING INITIAL POSITION...\n" << std::endl;

    std::cout << "\t\t -- Perifocal RF --\n" << std::endl;
    std::cout << "\t\tr0_p = [" << r_perifocal[0] << ", " << r_perifocal[1]
              << ", " << r_perifocal[2] << "]" << std::endl;
    std::cout << "\t\tv0_p = [" << v_perifocal[0] << ", " << v_perifocal[1]
              << ", " << v_perifocal[2] << "]\n"
              << std::endl;
    //
    std::cout << "\t\t -- ECI --\n" << std::endl;
    std::cout << "\t\tr0_ECI = [" << r_ECI[0] << ", " << r_ECI[1] << ", "
              << r_ECI[2] << "]" << std::endl;

    std::cout << "\t\tv0_ECI = [" << v_ECI[0] << ", " << v_ECI[1] << ", "
              << v_ECI[2] << "]\n"
              << std::endl;
        std::cout << "\t\t----------------------------------" << std::endl;   
   
    // Initialise
    double temp = 0.0;
    DU = 0.0;
    TU = 0.0;

    for (const auto &x : r_ECI) {
      temp += x * x;
    }
    DU = std::sqrt(temp);
    TU = std::sqrt((DU * DU * DU) / Earth_mu);

    // Initialise
    r_0 = std::vector<double>(3, 0.0);
    v_0 = std::vector<double>(3, 0.0);

    // Undimensionalise coordinates
    for (int i = 0; i < 3; i++) {
      r_0[i] = r_ECI[i] / DU;
      v_0[i] = v_ECI[i] / (DU / TU);
    }
    x = r_0;
    x.insert(x.end(), v_0.begin(), v_0.end());
};







std::vector<double> Orbit::ECI2LVLH(std::vector<double>&ECI){
    
    std::vector<double> ECI_r = {0.0,0.0,0.0};
    std::vector<double> ECI_v = {0.0,0.0,0.0};
    for(int i = 0; i < 3; i++){
        ECI_r[i] = ECI[i];
        ECI_v[i] = ECI[i+3];
    }

    double r_norm = sqrt(ECI_r[0] * ECI_r[0] + ECI_r[1] * ECI_r[1] + ECI_r[2] * ECI_r[2]);
    std::vector<double> LV = {0.0, 0.0, 0.0};
    
    for (int i = 0; i < 3; ++i){
        LV[i] = ECI_r[i]/r_norm;
    };

    
    std::vector<double> h(3,0.0);
    h = {(ECI_r[1]*ECI_v[2]-ECI_r[2]*ECI_v[1]), (ECI_r[0]*ECI_v[2]-ECI_r[2]*ECI_v[0]), (ECI_r[0]*ECI_v[1]-ECI_r[1]*ECI_v[0])};

    
    double h_norm = sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]);
    std::vector<double> orbit_normal = {0.0, 0.0, 0.0};
    for (int i = 0; i<3; ++i){
        orbit_normal[i] = h[i]/h_norm;
    };


    std::vector<double> LH(3,0.0); 
    LH = {(orbit_normal[1]*LV[2]-orbit_normal[2]*LV[1]), (orbit_normal[0]*LV[2]-orbit_normal[2]*LV[0]), (orbit_normal[0]*LV[1]-orbit_normal[1]*LV[0])};

    std::vector<std::vector<double>> ECILVLH = {
                {LV[0], LH[0], orbit_normal[0]},
                {LV[1], LH[1], orbit_normal[1]},
                {LV[2], LH[2], orbit_normal[2]},
                                                };
    std::vector<double> LVLH_x = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++){
            LVLH_x[i] += ECILVLH[i][j] * ECI_r[j];
            LVLH_x[i+3] += ECILVLH[i][j] * ECI_v[j];
        }
    }

    return LVLH_x;
}

std::vector<double> quaternion(std::vector<double> LVLH){

    std::vector<double> r = {LVLH[0], LVLH[1], LVLH[2]};
    double norm = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

    std::vector<double> LVLH_normalised = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; i++){
        LVLH_normalised[i] = r[i]/norm;
    };
    
    // pitch
    double phi = -asin(LVLH_normalised[0]);
    
    // roll
    double theta = atan2(LVLH_normalised[1], LVLH_normalised[2]);
    
    // yaw
    double psi = atan2(LVLH_normalised[1], LVLH_normalised[0]);
    

    // quaternion
    double qw = cos(phi/2.0);
    double qx = sin(phi/2.0) * cos(theta/2.0);
    double qy = sin(phi/2.0) * sin(theta/2.0) * cos(psi/2.0);
    double qz = sin(phi/2.0) * sin(theta/2.0) * sin(psi/2.0);
    
    std::vector<double> quaternion = {qw,qx,qy,qz};

    return quaternion;
}



/**
 * @brief a simplified approximation on GMST. 
 *
 */
double Orbit::GMST(double currentTime) {
    // You need to implement the GMST calculation here based on your specific needs.
    // This is a simplified example, and you may need a more accurate formula.
    // Here, we assume GMST increases linearly with time.
    const double gmstRate = 360.0 / 86164.100352; // Earth's rotation period (seconds)
    return fmod(280.46061837 + gmstRate * currentTime, 360.0);
}



/**
 * @brief Perform coordinate transformation from ECI to ECEF
 *
 * @param ECI:      Takes ECI vector
 * @param GMST_deg: Output from GMST function to construct rotational matrix
 */
std::vector<double> Orbit::EciToEcef(std::vector<double>& ECI, double GMST_deg){
    double GMST_rad = GMST_deg * (M_PI/180);

    std::vector<std::vector<double>> RotMat_Eci2Ecef = 
        {{cos(GMST_rad), -sin(GMST_rad), 0}, 
        {sin(GMST_rad), cos(GMST_rad), 0}, 
        {0, 0, 1}};
    
    std::vector<double> r_ECEF(3, 0.0);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            r_ECEF[i] += RotMat_Eci2Ecef[i][j] * ECI[j];
        }
    }
    return r_ECEF;
};



/** 
 * @brief Perform coordinate transformation from ECEF to Geodetic
 *
 * @param x: Takes ECEF vector
 * @param a: Semi major axis of Earth (approx.)
 * @param b: Semi minor acis of Earth (approx.)
 *
 */
std::vector<double> Orbit::EcefToGeo(std::vector<double>& ECEF){

    for (double& element : ECEF) {
        element*= DU; // Multiply each element with the result
    }
   
    const double a = 6378.137; // Semi-major axis in km
    const double b = 6356.7523; // Semi-minor axis in km
    //
    const double e2 = 1.0 - (b * b) / (a * a);

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
    
    // double r = sqrt(pow(ECEF[0],2.0) + pow(ECEF[1],2.0));
    // double E = (b*ECEF[2] - (a*a - b*b)) / (a*r);
    // double F = (b*ECEF[2] + (a*a - b*b)) / (a*r);
    // double P = (4.0/3.0) * (E*F + 1);
    // double Q = 2 * (E*E - F*F);
    // double D = pow(P,3.0) + pow(Q, 2.0);
    // 
    // double nu = 0.0;
    //
    // if (D < 0) {
    //     nu = 2*sqrt(P) * cos((1.0/3.0) * acos(Q*pow(-P,-0.5) / P));
    // }
    // else {
    //     nu = pow(pow(D,0.5) - Q, 1/3) - pow(pow(D,0.5) + Q , 1/3);
    // }
    //
    // double G = 0.5 * (sqrt(E*E + nu) + E);
    //
    // double temp_term = G*G + (F - nu*G) / (2*G - E);
    // double t = sqrt(temp_term) - G;
    //
    // double term = a * (1-t*t) / (2*b*t);
    
    return {longitude, latitude};
}



/**
 * @brief To contruct the Equation of motion of the satellite for integration accounting for gravity accerleration and oblateness accerleration
 *
 * @param x: state vector with component [r,v]
 *
 * @output x_dot: derivative of the state vector [v,a]
 */
std::vector<double> Orbit::EoM(std::vector<double> &x) {
    const double R_E = Earth_Radius / SMA;
    const double mu = 1.0;
    const double r_norm = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);

    // calculate and aggregate the acceleration from gravity and J2
    a_gravity = std::vector<double>(3, 0.0);
    a_oblateness = std::vector<double>(3, 0.0);

    const double common_term = -(3 * J2 * mu * R_E * R_E) / (2 * std::pow(r_norm, 5));
    const double z_squared_term = (5 * x[2] * x[2]) / (r_norm * r_norm);

    a_oblateness[0] = common_term * x[0] * (1 - z_squared_term);
    a_oblateness[1] = common_term * x[1] * (1 - z_squared_term);
    a_oblateness[2] = common_term * x[2] * (3 - z_squared_term);

    for (int i = 0; i < 3; ++i) {
      a_gravity[i] = -mu * x[i] / pow(r_norm, 3);
    }

    x_dot = std::vector<double>(6, 0.0);

    for (int i = 0; i < 3; ++i) {
      x_dot[i] = x[i + 3];
      x_dot[i + 3] = a_gravity[i] + a_oblateness[i];
    }

    return x_dot;
};



/**
 * @brief RungeKutta45 integration function
 *
 *
 */
void Orbit::RungeKutta45(double dt, double T, std::vector<double> &x) {

    std::vector<double> k1 = std::vector<double>(6, 0.0);
    std::vector<double> k2 = std::vector<double>(6, 0.0);
    std::vector<double> k3 = std::vector<double>(6, 0.0);
    std::vector<double> k4 = std::vector<double>(6, 0.0);
    std::vector<double> temp = std::vector<double>(6, 0.0);

    double timestep = T / dt;

    // Display Information on output text file
    std::cout << "\t\tWriting output to file..." << std::endl;
    std::ofstream vOut_ECI("ECI.txt", std::ios::out | std::ios::trunc);
 
    std::ofstream vOut_ECEF("ECEF.txt", std::ios::out | std::ios::trunc);
    std::ofstream vOut_GEO("GEO.txt", std::ios::out | std::ios::trunc);
    std::ofstream vOut_q("q.txt", std::ios::out | std::ios::trunc);
    

    for (int t = 0; t < timestep; t++) {

        k1 = Orbit::EoM(x);

        for (int i = 0; i < 6; i++)
            temp[i] = x[i] + (dt * k1[i]) / 2;
        
        k2 = Orbit::EoM(temp);

        for (int i = 0; i < 6; i++)
            temp[i] = x[i] + (dt * k2[i]) / 2;
        
        k3 = Orbit::EoM(temp);

        for (int i = 0; i < 6; i++)
            temp[i] = x[i] + dt * k3[i];
        k4 = Orbit::EoM(temp);
        
        for (int i = 0; i < 6; i++) {
            x[i] += (1.0 / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * dt;
        }


        // Checking that file opened successfully.
        if (vOut_ECI.is_open()) {
            vOut_ECI << x[0]*DU << " " << x[1]*DU << " " << x[2]*DU << " " << x[3] << " "
                << x[4] << " " << x[5] << std::endl;
        }
        else {
            std::cout << "Did not open file successfully!" << std::endl;
        }
        
        std::vector<double> LVLH(6,0.0), q(4,0.0), temp2(6,0.0);
        std::vector<double> temp1(3,0.0), ECEF(3,0.0);
        for(int i = 0; i < 3; i++){
            temp1[i] = x[i];
            temp2[i] = x[i];
            temp2[i+3] = x[i+3];
        }

        double current_time = Orbit::GMST(t*100); 
        ECEF = Orbit::EciToEcef(temp1, current_time);
        vOut_ECEF << ECEF[0]*DU << " " << ECEF[1]*DU << " " << ECEF[2]*DU << std::endl;



        LVLH = Orbit::ECI2LVLH(temp2);
        q = quaternion(LVLH);
        vOut_q << t << " " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << std::endl; 


        std::vector<double> GEO(2,0.0);
        GEO = Orbit::EcefToGeo(ECEF);
        vOut_GEO << GEO[0] << " " << GEO[1] << std::endl;      

        if (t % 628 == 0) {
            std::cout << "\t\tTime-step: " << std::setw(6) << t << "/" << timestep
                      << " (" << std::setw(2) << 100 * t / timestep << "%)"
                      << std::endl;
        }

   }
        
    // Finish writing and close file
    vOut_ECI.close();
    vOut_ECEF.close();
    vOut_GEO.close();
    vOut_q.close();
    std::cout << "\t\tFinished writing to file." << std::endl;
}



/**
 * @brief Perform integration with RungeKutta45 and print out the final position and velocity 
 *
 */
void Orbit::integrate() {

    // define target time and dt
    T = 11 * 2 * M_PI * sqrt(pow(SMA, 3.0) / Earth_mu) / TU;
    double dt = 0.01;

    // Initialise r_final and v_final vectors
    r_final = std::vector<double>(3, 0.0);
    v_final = std::vector<double>(3, 0.0);

    // Begin RK45 Integration
    Orbit::RungeKutta45(dt, T, x);

    // Dimensionalise r and v vectors
    for (int i = 0; i < 3; ++i) {
      r_final[i] = x[i] * DU;
      v_final[i] = x[i + 3] * (DU / TU);
    }

    std::cout << "\n\t\trf = [" << r_final[0] << ", " << r_final[1] << ", "
              << r_final[2] << "]" << std::endl;

    std::cout << "\t\tvf = [" << v_final[0] << ", " << v_final[1] << ", "
              << v_final[2] << "]\n"
              << std::endl;
    std::cout << "\t\t----------------------------------" << std::endl;
}



/**
 * @brief Perform 3x3 matrix multiplication
 *
 * @return std::vector<std::vector<double> >
 */
std::vector<std::vector<double>>
Orbit::matrixmultiply(std::vector<std::vector<double>> A,
                      std::vector<std::vector<double>> B) {
    std::vector<std::vector<double>> result(A.size(),
                                            std::vector<double>(B[0].size(), 0));
    if (A[0].size() == B.size()) {
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < B[0].size(); j++) {
            for (int k = 0; k < A[0].size(); k++)
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
};
//
// Orbit::~Orbit() {
//   r_perifocal.clear();
//   v_perifocal.clear();
//   r_ECI.clear();
//   v_ECI.clear();
//   r_0.clear();
//   v_0.clear();
//   x.clear();
//   x_dot.clear();
//   a_total.clear();
//   a_gravity.clear();
//   a_oblateness.clear();
//   r_final.clear();
//   v_final.clear();
//   RotMatrix_P2ECI.clear();
// };
