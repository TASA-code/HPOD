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
        const std::string Start_Date = "29-Oct-2023 07:55:48.000";
        const std::string End_Date = "30-Oct-2023 19:55:48.000";

        const double dt = 0.015625;
    
    public:

        static double SMA;
        static double e;
        static double i;
        static double RAAN;
        static double w;
        static double M;
        double h;
        double T;
        
        Vector6d ECI_state;
        Vector6d state;

        // Define the initial conditions from the uesr input
        void SetParameter(const double &arg_SMA, const double &arg_e, const double &arg_i,
                          const double &arg_M, const double &arg_w,
                          const double &arg_RAAN);
        
        // Equation of motion
        Vector6d EoM(Vector6d &x);

        void RungeKutta45(double dt, double T, Vector6d& x);

        void Propagate();

        // ~Orbit();
};
