#include <cmath>
#include <numeric>
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>

using namespace Eigen;

typedef Matrix<double,6,1> Vector6d;


/**
 * @class
 *
 */
class Orbit {

    private:
        const long double Earth_Radius = 6378.1;
        const long double Earth_mu = 398600.4415;
        const long double J2 = 0.00108263;
        
        const double dt = 0.015625;
    
    public:

        inline static std::string Start_Date;
        inline static std::string End_Date;
        inline static double SMA;
        inline static double e;
        inline static double i;
        inline static double RAAN;
        inline static double w;
        inline static double M;
        
        Vector6d state;

        // Define the initial conditions from the uesr input
        void SetParameter(const double &arg_SMA, const double &arg_e, const double &arg_i,
                          const double &arg_M, const double &arg_w,
                          const double &arg_RAAN, const std::string& arg_Start_Date, const std::string& arg_End_Date);
        
        // Equation of motion
        Vector6d f(const Vector6d &x);

        void RungeKutta45(const double& dt, const double& T, Vector6d& x);

        void Propagate();

        // ~Orbit();
};
