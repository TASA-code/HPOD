#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <ostream>
#include <vector>       
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>

#include "Orbit.h"
#include "Coordinate.h"



typedef Eigen::Matrix<double,6,1> Vector6d;


// static function initialisation
double Orbit::SMA   = 0.0;
double Orbit::e     = 0.0; 
double Orbit::i     = 0.0; 
double Orbit::RAAN  = 0.0;
double Orbit::w     = 0.0;
double Orbit::M     = 0.0;

double Orbit::TU = 0.0;
double Orbit::DU = 0.0;

/**
 * @brief Takes user input orbital element and set initial coordinates
 *
 * @param arg_SMA
 * @param arg_e
 * @param arg_i
 * @param arg_RAAN
 * @param arg_w
 * @param arg_M
 */
void Orbit::SetParameter(const double &arg_SMA, const double &arg_e,
                         const double &arg_i, const double &arg_RAAN,
                         const double &arg_w, const double &arg_M) {

    // Printing values of all the parameters to the terminal.
    std::cout << std::endl;
    std::cout << "\t\tList of Commands Line Inputs:\n" << std::endl;
    std::cout << "\t\t- SMA       :  " << arg_SMA << std::endl;
    std::cout << "\t\t- e         :  " << arg_e << std::endl;
    std::cout << "\t\t- i         :  " << arg_i << std::endl;
    std::cout << "\t\t- RAAN      :  " << arg_RAAN << std::endl;
    std::cout << "\t\t- w         :  " << arg_w << std::endl;
    std::cout << "\t\t- M         :  " << arg_M << std::endl;
    std::cout << std::endl;
    

    SMA = arg_SMA;
    e = arg_e;
    i = arg_i * M_PI / 180;
    RAAN = arg_RAAN * M_PI / 180;
    w = arg_w * M_PI / 180;
    M = arg_M * M_PI / 180;
    

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
    
    r_vector << cos(theta), sin(theta), 0;
    v_vector << -sin(theta), e+cos(theta), 0;

    P_r << (SMA*(1-e*e) / (1+e*cos(theta))) * r_vector;

    P_v << sqrt(Earth_mu/SMA) * v_vector;


    // calculate angular momentum h
    // h = sqrt(Earth_mu * SMA * (1 - e * e));
   
    // P_r << ( (h*h) / (Earth_mu * (1+e+cos(theta))) ) * (cos(theta) * i_e + sin(theta) * i_p);
    // P_v = (Earth_mu / h) * ( -sin(theta) * i_e + (e + cos(theta)) * i_p );

    Vector6d Perifocal;
    Perifocal << P_r, P_v;

    ECI_state = Coordinate::P2ECI(Perifocal); 

    Eigen::Vector3d ECI_r = ECI_state.head<3>();
    Eigen::Vector3d ECI_v = ECI_state.tail<3>();

    ECI_r << 918.8051813, 100.2960181, 6871.8919169;
    ECI_v << 0.884039, -7.534891, -0.0081953;
    

    // Print the resulting vectors r_ECI and v_ECI
    // format data output
    std::cout << "\t\t----------------------------------" << std::endl;
    std::cout << "\t\t...PRINTING INITIAL POSITION (ECI)...\n" << std::endl;
    std::cout << "\t\tr0 = [" << ECI_r[0] << ", " << ECI_r[1] << ", "
              << ECI_r[2] << "]" << std::endl;

    std::cout << "\t\tv0 = [" << ECI_v[0] << ", " << ECI_v[1] << ", "
              << ECI_v[2] << "]\n"
              << std::endl;
        std::cout << "\t\t----------------------------------" << std::endl;   
   

    // Initialise
    
    DU = ECI_r.norm();
    TU = sqrt((DU * DU * DU) / Earth_mu);

    // Initialise
    Eigen::Vector3d undimensional_r, undimensional_v;
    // Undimensionalise coordinates
    undimensional_r = ECI_r.array() / DU;
    undimensional_v = ECI_v.array() / (DU/TU);
    undimensional_state = (Vector6d(6) << undimensional_r, undimensional_v).finished();

}


/**
 * @brief To contruct the Equation of motion of the satellite for integration accounting for gravity accerleration and oblateness accerleration
 *
 * @param x: state vector with component [r,v]
 *
 * @output x_dot: derivative of the state vector [v,a]
 */
Vector6d Orbit::EoM(Vector6d& x) {

    Vector3d r = x.head<3>();
    Vector3d v = x.tail<3>();
    
    const double undimensional_mu = 1.0;
    const double r_norm = r.norm();


    // calculate and aggregate the acceleration from gravity and J2
    const double undimensional_RE = Earth_Radius / SMA;
    const double common_term = -(3 * J2 * undimensional_mu * undimensional_RE * undimensional_RE) / (2 * std::pow(r_norm, 5));
    const double z_squared_term = (5 * x[2] * x[2]) / (r_norm * r_norm);


    Eigen::Vector3d a_gravity, a_oblateness;
    a_oblateness[0] = common_term * x[0] * (1 - z_squared_term);
    a_oblateness[1] = common_term * x[1] * (1 - z_squared_term);
    a_oblateness[2] = common_term * x[2] * (3 - z_squared_term);

    a_gravity << -undimensional_mu*r / pow(r_norm,3);

    Vector6d x_dot;
    x_dot << v, a_gravity + a_oblateness;

    return x_dot;
};



/**
 * @brief RungeKutta45 integration function
 *
 *
 */
void Orbit::RungeKutta45(double dt, double T, Vector6d& x) {


    const int timestep = static_cast<int>(T/dt);
    // const int max_iteration = 2000000;

    int iteration = 0;
    double current_time = 0.0;


    // Display Information on output text file
    std::cout << "\t\tWriting output to file..." << std::endl;
    std::cout << std::endl;
    std::ofstream vOut_ECI("ECI.txt", std::ios::out | std::ios::trunc);
    std::ofstream vOut_ECEF("ECEF.txt", std::ios::out | std::ios::trunc);
    std::ofstream vOut_GEO("GEO.txt", std::ios::out | std::ios::trunc);
    std::ofstream vOut_q("q.txt", std::ios::out | std::ios::trunc);

    if (!vOut_ECI.is_open() || !vOut_ECEF.is_open() || !vOut_GEO.is_open() || !vOut_q.is_open()) {
        std::cout << "Failed to open output files!" << std::endl;
        return; // or handle the error as needed
    }
   


    // Initialise Vector
    Vector6d k1, k2, k3, k4, temp;

    // Perform Runge-Kutta algorithm
    while (current_time < T ) {

        k1 = Orbit::EoM(x);

        temp = x + (dt * k1) / 2.0;
        k2 = Orbit::EoM(temp);

        temp = x + (dt * k2) / 2.0;
        k3 = Orbit::EoM(temp);

        temp = x + dt * k3;
        k4 = Orbit::EoM(temp);

        Vector6d delta_x = (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4) * dt;

        // Check if the maximum absolute value in delta_x is less than the tolerance
        if (delta_x.cwiseAbs().maxCoeff() < 1e-15 ) {
            break;  // Exit the loop if the tolerance is satisfied
        }

        x += delta_x;
        current_time += dt;
        iteration++;

        

        // Write output to file (ECI, ECEF, GEO)
        Eigen::VectorXd Output_ECI(7);
        Output_ECI << current_time*TU, x.head<3>() * DU, x.segment<3>(3) * (DU / TU);
        vOut_ECI << Output_ECI.transpose() << std::endl;
        
        // Transform ECI to ECEF and output to file 
        Vector6d ECEF;  
        ECEF = Coordinate::ECI2ECEF(x, current_time*TU);
        vOut_ECEF << ECEF[0]*DU << " " << ECEF[1]*DU << " " 
                  << ECEF[2]*DU << std::endl;

        // Transform ECEF to GEO and output to file
        Vector2d GEO; 
        GEO = Coordinate::ECEF2GEO(ECEF);
        vOut_GEO << GEO[0] << " " << GEO[1] << std::endl;      
     


        // Progress Bar Display
        double progress = static_cast<double>(iteration) / timestep;
        const int barWidth = 65;
        int pos = static_cast<int>(barWidth * progress);

        std::cout << "[";
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos)
                std::cout << "#";
            else if (i == pos)
                std::cout << "#";
            else
                std::cout << " ";
        }
        
        std::cout << "] " << std::fixed << std::setprecision(3)
                        << progress * 100.0 << "%\r";
        std::cout.flush();


    // for (int t = 0; t < timestep; t++) {
    // 
    //     k1 << Orbit::EoM(x);
    //
    //     temp << x + (dt * k1) / 2.0;
    //
    //     k2 << Orbit::EoM(temp);
    //     
    //     temp << x + (dt * k2) / 2.0;
    //     
    //     k3 << Orbit::EoM(temp);
    //     
    //     temp << x + dt * k3;
    //
    //     k4 << Orbit::EoM(temp);
    //     
    //     x += (1.0/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4) * dt;
    //     

       //  Eigen::VectorXd Output_ECI(7);
       //  Output_ECI << t, x.head<3>() * DU, x.segment<3>(3) * (DU / TU);
       //  vOut_ECI << Output_ECI.transpose() << std::endl;
       //  
       //  
       //  Vector6d ECEF;  
       //  // double current_time = Coordinate::GMST(t); // t*100 
       //  ECEF = Coordinate::ECI2ECEF(x, t);
       //  vOut_ECEF << ECEF[0]*DU << " " << ECEF[1]*DU << " " << ECEF[2]*DU << std::endl;
       //
       // 
       //  Vector2d GEO; 
       //  GEO = Coordinate::ECEF2GEO(ECEF);
       //  vOut_GEO << GEO[0] << " " << GEO[1] << std::endl;      
       //
       //  if (t % 628 == 0) {
       //      std::cout << "\t\tTime-step: " << std::setw(6) << t << "/" << timestep
       //                << " (" << std::setw(2) << 100 * t / timestep << "%)"
       //                << std::endl;
       //  }
   }
        
    // Finish writing and close file
    vOut_ECI.close();
    vOut_ECEF.close();
    vOut_GEO.close();
    vOut_q.close();
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "\t\tFinished writing to file." << std::endl;
}



/**
 * @brief Perform integration with RungeKutta45 and print out the final position and velocity 
 *
 */
void Orbit::integrate() {

    // define target time and dt
    const double T = 13 * 2 * M_PI * sqrt(pow(SMA, 3.0) / Earth_mu) / TU; //75,478.687
    const double dt = 0.01/TU;

    // Begin RK45 Integration
    Orbit::RungeKutta45(dt, T, undimensional_state);
    
    // Initialise final state vector
    // Dimensionalise r and v vectors
    Vector3d Final_r = undimensional_state.head<3>() * DU;
    Vector3d Final_v = undimensional_state.tail<3>() * (DU/ TU);

    std::cout << "\n\t\trf = [" << Final_r[0] << ", " << Final_r[1] << ", "
              << Final_r[2] << "]" << std::endl;

    std::cout << "\t\tvf = [" << Final_v[0] << ", " << Final_v[1] << ", "
              << Final_v[2] << "]\n"
              << std::endl;
    std::cout << "\t\t----------------------------------" << std::endl;
}

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
