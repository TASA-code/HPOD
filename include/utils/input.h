#ifndef _READINPUT_H
#define _READINPUT_H

#include <string>

struct SatelliteData {
    double SMA = 0.0;
    double e = 0.0;
    double i = 0.0;
    double M = 0.0;
    double w = 0.0;
    double RAAN = 0.0;
    std::string Start_Date = " ";
    std::string End_Date = " ";
};

SatelliteData parseInputFile(const std::string& filePath);

#endif // MYSTRUCT_H