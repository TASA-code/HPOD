#include <iomanip>
#include <iostream>
#include <fstream>

#include "propagator.h"
#include "accel.h"
#include "coordinate.h"
#include "dm_time.h"
#include "rk45.h"


#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>
using namespace Eigen;
typedef Eigen::Matrix<double,6,1> Vector6d;


// Progress bar function
void Progress_Bar(int& iteration, const int& timestep){
    double progress = static_cast<double>(iteration) / timestep;
    const int barWidth = 35;
    int pos = static_cast<int>(barWidth * progress);
    
    // ANSI escape code for green
    const std::string reset = "\033[0m";
    const std::string green = "\033[32m";

    std::cout << "\t\t[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos)
            std::cout << green << "#" << reset;
        else if (i == pos)
            std::cout << green << "#" << reset;
        else
            std::cout << " ";
    }
    
    std::cout << "] " << std::fixed << std::setprecision(3)
                    << progress * 100.0 << "%\r";
    std::cout.flush();
}




void WriteOutputToFile(const std::string& date, const Vector6d& x, const Vector6d& ECEF, const Vector2d& GEO1, const Vector2d& GEO2) {

    std::ofstream outputFile("output.txt", std::ios::out | std::ios::app);

    if (!outputFile.is_open()) {
        std::cout << "Failed to open output files!" << std::endl;
        return; // or handle the error as needed
    }

    // Set a fixed width for each field (adjust the widths as needed)
    const int fieldWidth = 5;

    // Write each field with a fixed width and a space as a separator
    outputFile << std::left << date;
    outputFile << std::setw(fieldWidth) << std::left << x.transpose();
    outputFile << std::setw(fieldWidth) << std::left << ECEF.transpose();
    outputFile << std::setw(fieldWidth) << std::left << GEO1.transpose();
    outputFile << std::setw(fieldWidth) << std::left << GEO2.transpose() << std::endl;  // Append a newline


    // Ensure the data is immediately written to the file
    outputFile.flush();

    outputFile.close();
}







void RungeKutta45(const double& T, const double& dt, const int& outputFrequency, Vector6d& x)
{

    const int timestep      = static_cast<int>(T/dt);
    int iteration           = 0;
    double current_time     = 0.0;
    
    std::string date;
    Vector6d k1, k2, k3, k4, delta_x; 
    Vector6d ECEF;
    Vector2d GEO1, GEO2;


    // Display Information on output text file
    std::cout << "\t\t...ORBIT PROPAGATING..." << std::endl << std::endl;
    
    // Perform Runge-Kutta algorithm
    while (current_time < T ) {

        k1 = f(x, current_time);

        k2 = f((x + (dt * k1) / 2.0), current_time);
        
        k3 = f((x + (dt * k2) / 2.0), current_time);

        k4 = f(x + dt * k3, current_time);

        delta_x.noalias() = (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4) * dt;


        // Write output to file every 'outputFrequency' iterations
        if (iteration % outputFrequency == 0) {
            // Write output to file (ECI, ECEF, GEO)
            date = Time2Date(Propagator::Start_Date, current_time);

            // Transform ECI to ECEF and output to file 
            ECEF = ECI2ECEF(x, current_time);  

            // Transform ECEF to GEO and output to file
            GEO1 = ECEF2GEO(x);
            GEO2 = ECEF2GEO(ECEF); 
            
            WriteOutputToFile(date, x, ECEF, GEO1, GEO2);
        };
        // Progress Bar Display
        Progress_Bar(iteration, timestep);
        

        x += delta_x;
        current_time += dt;
        iteration++;
   }
        
    // Finish writing and close file
    std::cout << "\n" << std::endl;
    std::cout << "\t\tPropagation Complete." << std::endl;
}


