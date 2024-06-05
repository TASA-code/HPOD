#include <iomanip>
#include <iostream>
#include <cmath>
#include <ctime>
#include <chrono>
#include <sstream>
#include <map>


#include "time/time.h"
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>
using namespace Eigen;




Eigen::VectorXi Extract_Date(const std::string& Date) {
std::istringstream iss(Date);

    int day, month, year, hour, minute, second, milliseconds;
    char dot; // Variable to store the dot separator
    std::string monthstr;

    iss >> day;
    iss.ignore(); // ignore the hyphen
    std::getline(iss, monthstr, '-'); // read until the next hyphen
    iss >> year >> hour >> dot >> minute >> dot >> second >> dot >> milliseconds;

    // Map month names to numerical values
    std::map<std::string, int> monthMap = {
        {"Jan", 1}, {"Feb", 2}, {"Mar", 3}, {"Apr", 4}, {"May", 5}, {"Jun", 6},
        {"Jul", 7}, {"Aug", 8}, {"Sep", 9}, {"Oct", 10}, {"Nov", 11}, {"Dec", 12}
    };

    // Convert month string to numerical value
    auto it = monthMap.find(monthstr);
    if (it != monthMap.end()) {
        month = it->second;
    } else {
        // Handle invalid month case
        std::cerr << "Invalid month: " << monthstr << std::endl;
        // You might want to handle this case differently, depending on your requirements
    }

    Eigen::VectorXi intDate(7); // Include milliseconds
    intDate << day, month, year, hour, minute, second, milliseconds;

    return intDate;
}



std::string Time2Date(const std::string& initialTime, double secondsToAdd) {
    // Parse the initial time string
    struct std::tm tm = {};
    std::istringstream ss(initialTime);

    // Modify the time format to include milliseconds during parsing
    ss >> std::get_time(&tm, "%d-%b-%Y %H:%M:%S");

    // Convert to std::chrono::system_clock time_point
    auto tp = std::chrono::system_clock::from_time_t(std::mktime(&tm));

    // Add seconds (including fractional part)
    tp += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration<double>(secondsToAdd));

    // Convert back to std::tm
    std::time_t tt = std::chrono::system_clock::to_time_t(tp);
    std::tm new_tm = *std::localtime(&tt);

    // Format the new time as a string
    std::ostringstream result;
    result << std::put_time(&new_tm, "%d-%b-%Y %H:%M:%S");

    // Include milliseconds
    auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(tp.time_since_epoch()).count() % 1000;
    result << "." << std::setw(3) << std::setfill('0') << milliseconds;

    return result.str();
}




double Duration(const std::string& startDate, const std::string& endDate) {
    // Parse the start date string
    struct std::tm tmStart = {};
    std::istringstream ssStart(startDate);
    ssStart >> std::get_time(&tmStart, "%d-%b-%Y %H:%M:%S");

    // Parse the end date string
    struct std::tm tmEnd = {};
    std::istringstream ssEnd(endDate);
    ssEnd >> std::get_time(&tmEnd, "%d-%b-%Y %H:%M:%S");

    // Convert to std::chrono::system_clock time_point
    auto tpStart = std::chrono::system_clock::from_time_t(std::mktime(&tmStart));
    auto tpEnd = std::chrono::system_clock::from_time_t(std::mktime(&tmEnd));

    // Calculate the time difference in seconds
    std::chrono::duration<double> duration = tpEnd - tpStart;
    return duration.count();
}