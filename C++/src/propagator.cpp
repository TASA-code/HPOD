#include <iomanip>
#include <iostream>
#include <ostream>
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>
#include <chrono>

#include "propagator.h"
#include "coordinate/coordinate.h"
#include "time/time.h"
#include "accel/accel.h"
#include "utils/rk45.h"



typedef Eigen::Matrix<double,6,1> Vector6d;




/**
 * @brief Takes user input orbital element and set initial coordinates
 *
 */
void Propagator::Initialise(const double &arg_SMA, const double &arg_e,
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
    state = P2ECI(Perifocal); 

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
 * @brief Perform integration with RungeKutta45 and print out the final position and velocity 
 *
 */
void Propagator::Propagate() {

    // define target time and dt
    const double T = Duration(Start_Date,End_Date);

    // Begin RK45 Integration
    auto start = std::chrono::high_resolution_clock::now();
    RungeKutta45(T, dt, output_frequency, state);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);

    std::cout << "\n\t\tPropagation time: " << duration.count() << " seconds" << std::endl;

    
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

