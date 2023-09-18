#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

/**
 * @class
 *
 */
class Orbit {

private:
  const double Earth_Radius = 6378;
  const double Earth_mu = 398600;
  const double J2 = 0.00108263;

public:
  double SMA;
  double e;
  double i;
  double RAAN;
  double w;
  double theta;
  double h;
  double T;

  double DU;
  double TU;

  std::vector<double> r_p;
  std::vector<double> v_p;
  std::vector<double> r_ECI;
  std::vector<double> v_ECI;
  std::vector<double> r_0;
  std::vector<double> v_0;

  std::vector<double> x;
  std::vector<double> x_dot;
  std::vector<double> a;
  std::vector<double> ag;
  std::vector<double> ad;

  std::vector<double> x_next;

  std::vector<double> r_final;
  std::vector<double> v_final;

  std::vector<std::vector<double>> A;

  // Define the initial conditions from the uesr input
  void SetParameter(const double &arg_SMA, const double &e, const double &arg_i,
                    const double &arg_RAAN, const double &arg_w,
                    const double &arg_theta);

  // Coordinate transformation to ECI and undimensionalise
  void P2ECI();

  std::vector<double> EoM(std::vector<double> x);

  void integrate();

  std::vector<std::vector<double>>
  matrixmultiply(std::vector<std::vector<double>> A,
                 std::vector<std::vector<double>> B);

  ~Orbit();
};
