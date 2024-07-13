#include "propagator.h"
#include "coordinate.h"
#include "dm_time.h"
#include "input.h"

#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>



int main(int argc, char *argv[]) {

    SatelliteData data = parseInputFile("../../input.txt");

    Propagator prop;


    prop.Initialise(data.SMA, data.e, data.i, data.M, data.w, data.RAAN, 
                    data.Start_Date, data.End_Date, data.step_time, data.sample_rate);
    
    prop.Propagate();

    
    // orbit_prop.~Orbit();
}
