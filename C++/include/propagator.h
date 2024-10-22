
#include <numeric>
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>

using namespace Eigen;

typedef Matrix<double,6,1> Vector6d;


/**
 * @class
 *
 */
class Propagator {

    public:

        inline static double Earth_Radius = 6378.1363e3;
        inline static double Earth_mu = 398600.4415e+9;

        inline static std::string Start_Date;
        inline static std::string End_Date;
        inline static double SMA;
        inline static double e;
        inline static double i;
        inline static double RAAN;
        inline static double w;
        inline static double M;
        inline static double step_time;
        inline static int sample_rate;
        
        Vector6d state;

        // Define the initial conditions from the uesr input
        void Initialise(const double &arg_SMA, const double &arg_e, const double &arg_i, 
                            const double &arg_M, const double &arg_w, const double &arg_RAAN, const std::string& arg_Start_Date, 
                            const std::string& arg_End_Date, const double& arg_step_time, const int& arg_sample_rate);
        

        void Propagate();

        // ~Orbit();
};
