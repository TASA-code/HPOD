#include <iomanip>
#include <iostream>
#include <fstream>

#include "propagator.h"
#include "accel/accel.h"
#include "coordinate/coordinate.h"
#include "time/time.h"
#include "utils/rk45.h"


#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>
using namespace Eigen;
typedef Eigen::Matrix<double,6,1> Vector6d;


// Progress bar function
void Progress_Bar(int& iteration, const int& timestep){
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



void RungeKutta45(const double& T, const double& dt, Vector6d& x)
{

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
        std::string date = Time2Date(Propagator::Start_Date, current_time);

        vOut_ECI << date << " " << Output_ECI.transpose() << std::endl;
        
        // Transform ECI to ECEF and output to file 
        Vector6d ECEF;  
        ECEF = ECI2ECEF(x, current_time);
        vOut_ECEF << ECEF.transpose() << std::endl;      

        // Transform ECEF to GEO and output to file
        Vector2d GEO1; 
        GEO1 = ECEF2GEO(x);
        vOut_GEO1 << GEO1.transpose() << std::endl;    

        Vector2d GEO2;
        GEO2 = ECEF2GEO(ECEF); 
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


