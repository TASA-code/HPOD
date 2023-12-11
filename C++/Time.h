#include <iomanip>
#include <iostream>
#include <cmath>
#include <ctime>
#include <chrono>


class Time {

public:    
    static std::string Time2Date(const std::string& initialTime, double secondsToAdd);
    static double Duration(const std::string& startDate, const std::string& endDate);
};
