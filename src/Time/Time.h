#ifndef  _time_H
#define  _time_H

#include <iostream>
#include <chrono>

#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>


extern Eigen::VectorXi Extract_Date(const std::string& Date);
    
extern std::string Time2Date(const std::string& initialTime, double secondsToAdd);
    
extern double Duration(const std::string& startDate, const std::string& endDate);


#endif