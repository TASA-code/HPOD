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


    double OE[] = {SMA, e, i, M, w, RAAN};

    state = OE2ECI(OE);
    
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
    auto duration = std::chrono::duration<double>(stop - start);

    std::cout << "\n\t\tPropagation time: " << std::fixed << std::setprecision(3) << duration.count() << " seconds" << std::endl;

    
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

