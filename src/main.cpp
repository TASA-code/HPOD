#include "Orbit/Orbit.h"
#include "Coordinate/Coordinate.h"
#include "Time/Time.h"
#include <iostream>


int main(int argc, char *argv[]) {
  
    const double SMA = 6982;
    const double e = 0.0010293460;
    const double i = 97.97357;
    const double M = 229.9764;
    const double w = 98.03756;
    const double RAAN = 340.3470;

    Orbit orbit_prop;

    orbit_prop.SetParameter(SMA, e, i, M, w, RAAN);

    orbit_prop.Propagate();

    // orbit_prop.~Orbit();
}
 