#include <iomanip>
#include <iostream>
#include <cmath>
#include <ctime>
#include <chrono>
#include <ostream>
#include <fstream>
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>

#include "Orbit.h"
#include "../Coordinate/Coordinate.h"
#include "../Time/Time.h"
#include "../Potential/accel.h"
// #include "../Potential/Potential.h"


typedef Eigen::Matrix<double,6,1> Vector6d;

namespace {

    // Progress bar function
    static void Progress_Bar(int& iteration, const int& timestep){
        double progress = static_cast<double>(iteration) / timestep;
        const int barWidth = 55;
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
    }
}




/**
 * @brief Takes user input orbital element and set initial coordinates
 *
 */
void Orbit::SetParameter(const double &arg_SMA, const double &arg_e,
                         const double &arg_i, const double &arg_M,
                         const double &arg_w, const double &arg_RAAN, 
                         const std::string& arg_Start_Date, 
                         const std::string& arg_End_Date) {

    // Printing values of all the parameters to the terminal.
    std::cout << std::endl;
    std::cout << "\t\tList of Commands Line Inputs:\n" << std::endl;
    std::cout << "\t\t- SMA       :  " << arg_SMA << std::endl;
    std::cout << "\t\t- e         :  " << arg_e << std::endl;
    std::cout << "\t\t- i         :  " << arg_i << std::endl;
    std::cout << "\t\t- M         :  " << arg_M << std::endl;
    std::cout << "\t\t- w         :  " << arg_w << std::endl;
    std::cout << "\t\t- RAAN      :  " << arg_RAAN << std::endl;
    std::cout << std::endl;
    std::cout << "\t\tStart Time:  " << arg_Start_Date << std::endl;
    std::cout << "\t\tEnd Time  :  " << arg_End_Date << std::endl;
    std::cout << std::endl;

    Start_Date = arg_Start_Date;
    End_Date = arg_End_Date;
    
    SMA = arg_SMA;
    e = arg_e;
    i = arg_i * M_PI / 180;
    M = arg_M * M_PI / 180;
    w = arg_w * M_PI / 180;
    RAAN = arg_RAAN * M_PI / 180;
    
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


    Vector6d Perifocal;
    Perifocal << P_r, P_v;
    state = Coordinate::P2ECI(Perifocal); 

    // ECI_r << 5748.272127, -1506.348412, -3674.401874;
    // ECI_v << 3.472888, -2.183142, 6.339326;

    state << 5748.272127e3, -1506.348412e3, -3674.401874e3, 
             3.472888e3, -2.183142e3, 6.339326e3;
    

    // Print the resulting vectors r_ECI and v_ECI
    // format data output
    std::cout << "\t\t----------------------------------" << std::endl;
    std::cout << "\t\t...INITIAL POSITION (ECI)...\n" << std::endl;
    std::cout << "\t\tr0 = [" << state[0] << ", " << state[1] << ", "
              << state[2] << "]" << std::endl;

    std::cout << "\t\tv0 = [" << state[3] << ", " << state[4] << ", "
              << state[5] << "]\n"
              << std::endl;
    std::cout << "\t\t----------------------------------" << std::endl;   
}




/**
 * @brief To contruct the Equation of motion of the satellite for integration accounting for gravity accerleration and oblateness accerleration
 *
 */
Vector6d Orbit::f(const Vector6d& x, double t) {

    Vector3d r_GCRF, v_GCRF;
    r_GCRF = x.head<3>();
    v_GCRF = x.tail<3>();

    // Vector3d r_ITRF,a_ITRF;
    // Vector6d x_ITRF;
    // x_ITRF = Coordinate::ECI2ECEF(x,t);
    // r_ITRF = x_ITRF.head<3>();
    // Vector2d GEO;
    // GEO = Coordinate::ECEF2GEO(x_ITRF);
    // double phi, lambda;
    // lambda = GEO[0];
    // phi = GEO[1];
    // a_ITRF = Potential::a(r_ITRF, phi, lambda);
    
    Vector3d a_GCRF;
    Vector6d x_dot;

    // a_GCRF = Coordinate::ECEF2ECI(a_ITRF,t);

    a_GCRF = AccelMod(x, Earth_mu, 10, 10, Earth_Radius,t);
    x_dot << v_GCRF, a_GCRF;
    return x_dot; 
};



/**
 * @brief RungeKutta45 integration function
 *
 */
void Orbit::RungeKutta45(const double& T, Vector6d& x) {


    const int timestep = static_cast<int>(T/dt);
    int iteration = 0;
    double current_time = 0.0;


    // Display Information on output text file
    std::cout << "\t\t...ORBIT PROPAGATING..." << std::endl << std::endl;
    std::ofstream vOut_ECI("ECI.txt", std::ios::out | std::ios::trunc), 
                  vOut_ECEF("ECEF.txt", std::ios::out | std::ios::trunc), 
                  vOut_GEO1("ECI_GROUNDTRACK.txt", std::ios::out | std::ios::trunc),
                  vOut_GEO2("GROUNDTRACK.txt", std::ios::out | std::ios::trunc);
    vOut_ECI << std::fixed << std::setprecision(6);


    if (!vOut_ECI.is_open() || !vOut_ECEF.is_open() || !vOut_GEO1.is_open() || !vOut_GEO2.is_open()) {
        std::cout << "Failed to open output files!" << std::endl;
        return; // or handle the error as needed
    }
   
    // Initialise Vectorcd
    Vector6d k1, k2, k3, k4; 

    // Perform Runge-Kutta algorithm
    while (current_time < T ) {

        k1 = f(x, current_time);

        k2 = f((x + (dt * k1) / 2.0), current_time);
        
        k3 = f((x + (dt * k2) / 2.0), current_time);

        k4 = f(x + dt * k3, current_time);

        Vector6d delta_x = (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4) * dt;

        // Check if the maximum absolute value in delta_x is less than the tolerance
        if (delta_x.cwiseAbs().maxCoeff() < 1e-15 ) {
            break;  // Exit the loop if the tolerance is satisfied
        }


        // Write output to file (ECI, ECEF, GEO)
        Eigen::VectorXd Output_ECI(6);
        Output_ECI << x.head<3>() , x.segment<3>(3);
        std::string date = Time::Time2Date(Start_Date,current_time);

        vOut_ECI << date << " " << Output_ECI.transpose() << std::endl;
        
        // Transform ECI to ECEF and output to file 
        Vector6d ECEF;  
        ECEF = Coordinate::ECI2ECEF(x, current_time);
        vOut_ECEF << ECEF.transpose() << std::endl;      

        // Transform ECEF to GEO and output to file
        Vector2d GEO1; 
        GEO1 = Coordinate::ECEF2GEO(x);
        vOut_GEO1 << GEO1.transpose() << std::endl;    

        Vector2d GEO2;
        GEO2 = Coordinate::ECEF2GEO(ECEF); 
        vOut_GEO2 << GEO2.transpose() << std::endl;    
    
        // Progress Bar Display
        Progress_Bar(iteration, timestep);
        

        x += delta_x;
        current_time += dt;
        iteration++;
   }
        
    // Finish writing and close file
    vOut_ECI.close();
    vOut_ECEF.close();
    vOut_GEO1.close();
    vOut_GEO2.close();
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "\t\tFinished writing to file." << std::endl;
}



/**
 * @brief Perform integration with RungeKutta45 and print out the final position and velocity 
 *
 */
void Orbit::Propagate() {

    // define target time and dt
    const double T = Time::Duration(Start_Date,End_Date);

    // Begin RK45 Integration
    Orbit::RungeKutta45(T, state);
    
    // Initialise final state vector
    Vector3d Final_r = state.head<3>() ;
    Vector3d Final_v = state.tail<3>() ;

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
