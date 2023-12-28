#include "Orbit/Orbit.h"
#include "Coordinate/coordinate.h"
#include "Time/time.h"
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>


#include <iostream>
#include <fstream>
#include <regex>


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


SatelliteData parseInputFile(const std::string& filePath) {
    std::ifstream inputFile(filePath);

    if (!inputFile.is_open()) {
        std::cerr << "ERROR OPENING FILE" << std::endl;
        return SatelliteData(); // Return a default-constructed object
    }

    SatelliteData satelliteData;

    std::regex regexPatterns[] = {
        std::regex("SMA\\s*\\(km\\)\\s*:\\s*(\\d+)"),
        std::regex("e\\s*\\(ND\\)\\s*:\\s*([\\d.]+)"),
        std::regex("i\\s*\\(deg\\)\\s*:\\s*([\\d.]+)"),
        std::regex("M\\s*\\(deg\\)\\s*:\\s*([\\d.]+)"),
        std::regex("w\\s*\\(deg\\)\\s*:\\s*([\\d.]+)"),
        std::regex("RAAN\\s*\\(deg\\)\\s*:\\s*([\\d.]+)"),
        std::regex("Start_Date\\s*:\\s*([^\n]+)"),
        std::regex("End_Date\\s*:\\s*([^\n]+)")
    };

    std::smatch match;
    std::string line;

    for (std::size_t i = 0; i < sizeof(regexPatterns) / sizeof(regexPatterns[0]); ++i) {
        while (std::getline(inputFile, line)) {
            if (std::regex_search(line, match, regexPatterns[i])) {
                switch (i) {
                    case 0: satelliteData.SMA = std::stod(match[1].str()); break;
                    case 1: satelliteData.e = std::stod(match[1].str()); break;
                    case 2: satelliteData.i = std::stod(match[1].str()); break;
                    case 3: satelliteData.M = std::stod(match[1].str()); break;
                    case 4: satelliteData.w = std::stod(match[1].str()); break;
                    case 5: satelliteData.RAAN = std::stod(match[1].str()); break;
                    case 6: satelliteData.Start_Date = match[1].str(); break;
                    case 7: satelliteData.End_Date = match[1].str(); break;
                }
            }
        }
        inputFile.clear(); // Clear the end-of-file flag
        inputFile.seekg(0, std::ios::beg); // Reset file position to the beginning
    }

    inputFile.close();
    return satelliteData;
}



int main(int argc, char *argv[]) {

    SatelliteData data = parseInputFile("../../input.txt");

    Orbit orbit_prop;


    orbit_prop.SetParameter(data.SMA, data.e, data.i, data.M, data.w, data.RAAN, data.Start_Date, data.End_Date);

    orbit_prop.Propagate();

    // orbit_prop.~Orbit();
}
