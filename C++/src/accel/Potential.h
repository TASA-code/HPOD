#include <algorithm>
#include <cmath>
#include <numeric>
#include <map>
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>


using namespace Eigen;

typedef Matrix<double,6,1> Vector6d;

/**
* @class 
*
*/
class Potential {

    private:
        
        inline static double R = Orbit::Earth_Radius; 
        inline static double mu = Orbit::Earth_mu;
        inline static int l_max = 7;

    public:
    
        static double Legendre_poly(int l, int m, double gamma);

        static double dU_dr(Vector3d& r, double phi, double lambda);

        static double dU_dphi(Vector3d& r, double phi, double lambda); 

        static double dU_dlambda(Vector3d& r, double phi, double lambda);

        static Vector3d a(Vector3d& r, double phi, double lambda);

}; 
