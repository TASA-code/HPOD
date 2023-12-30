#include "propagator.h"
#include "coordinate/coordinate.h"
#include "time/time.h"
#include "utils/input.h"
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>



int main(int argc, char *argv[]) {

    SatelliteData data = parseInputFile("../../input.txt");

    Propagator prop;


    prop.SetParameter(data.SMA, data.e, data.i, data.M, data.w, data.RAAN, data.Start_Date, data.End_Date);

    prop.Propagate();

    // orbit_prop.~Orbit();
}
