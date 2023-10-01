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

  std::vector<double> r_perifocal;
  std::vector<double> v_perifocal;
  std::vector<double> r_ECI;
  std::vector<double> v_ECI;
  std::vector<double> r_0;
  std::vector<double> v_0;

  std::vector<double> x;
  std::vector<double> x_dot;
  std::vector<double> a_total;
  std::vector<double> a_gravity;
  std::vector<double> a_oblateness;

  std::vector<double> r_final;
  std::vector<double> v_final;

  std::vector<std::vector<double>> RotMatrix_P2ECI;

  // Define the initial conditions from the uesr input
  void SetParameter(const double &arg_SMA, const double &e, const double &arg_i,
                    const double &arg_RAAN, const double &arg_w,
                    const double &arg_theta);

  // Coordinate transformation to ECI and undimensionalise
  void P2ECI();

  // Coordinate transformation from ECI to ECEF
  std::vector<double> EciToEcef(std::vector<double> &x);

  std::vector<double> EoM(std::vector<double> &x);

  void RungeKutta45(double dt, double T, std::vector<double> &x);

  void integrate();

  std::vector<std::vector<double>>
  matrixmultiply(std::vector<std::vector<double>> A,
                 std::vector<std::vector<double>> B);
};
