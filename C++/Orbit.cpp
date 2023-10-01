#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include "Orbit.h"

/**
 * @brief Takes user input orbital element and set initial coordinates
 *
 * @param arg_SMA
 * @param e
 * @param arg_i
 * @param arg_RAAN
 * @param arg_w
 * @param arg_theta
 */
void Orbit::SetParameter(const double &arg_SMA, const double &arg_e,
                         const double &arg_i, const double &arg_RAAN,
                         const double &arg_w, const double &arg_theta) {

  SMA = arg_SMA;
  e = arg_e;
  i = arg_i * M_PI / 180;
  RAAN = arg_RAAN * M_PI / 180;
  w = arg_w * M_PI / 180;
  theta = arg_theta * M_PI / 180;

  // initialise vector
  r_perifocal = {0.0, 0.0, 0.0};
  v_perifocal = {0.0, 0.0, 0.0};
  std::vector<double> i_e = {1.0, 0.0, 0.0};
  std::vector<double> i_p = {0.0, 1.0, 0.0};

  // calculate angular momentum h
  h = std::sqrt(Earth_mu * SMA * (1 - e * e));

  for (int i = 0; i < 3; ++i) {
    r_perifocal[i] = ((h * h) / (Earth_mu * (1 + e * std::cos(theta)))) *
                     (std::cos(theta) * i_e[i] + std::sin(theta) * i_p[i]);
    v_perifocal[i] = (Earth_mu / h) * (-std::sin(theta) * i_e[i] +
                                       (e + std::cos(theta)) * i_p[i]);
  }

  // Printing values of all the parameters to the terminal.
  std::cout << std::endl;
  std::cout << "\t\tList of Commands Line Inputs:\n" << std::endl;
  std::cout << "\t\t- SMA       :  " << SMA << std::endl;
  std::cout << "\t\t- e         :  " << e << std::endl;
  std::cout << "\t\t- i         :  " << i << std::endl;
  std::cout << "\t\t- RAAN      :  " << RAAN << std::endl;
  std::cout << "\t\t- w         :  " << w << std::endl;
  std::cout << "\t\t- theta     :  " << theta << std::endl;
  std::cout << std::endl;
};

/**
 * @brief Perform coordinate transform from perifocal frame to ECI
 *
 */
void Orbit::P2ECI() {

  // Define transformation matrices R_RAAN, R_i, and R_w
  std::vector<std::vector<double>> R_RAAN = {
      {cos(RAAN), sin(RAAN), 0}, {-sin(RAAN), cos(RAAN), 0}, {0, 0, 1}};

  std::vector<std::vector<double>> R_i = {
      {1, 0, 0}, {0, cos(i), sin(i)}, {0, -sin(i), cos(i)}};

  std::vector<std::vector<double>> R_w = {
      {cos(w), sin(w), 0}, {-sin(w), cos(w), 0}, {0, 0, 1}};

  // Calculate the transformation matrix A = R_RAAN * R_i * R_w
  // Initialise matrix A
  RotMatrix_P2ECI =
      std::vector<std::vector<double>>(3, std::vector<double>(3, 0.0));
  RotMatrix_P2ECI = matrixmultiply(matrixmultiply(R_w, R_i), R_RAAN);

  std::vector<std::vector<double>> R_transpose(3, std::vector<double>(3, 0.0));
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      R_transpose[i][j] = RotMatrix_P2ECI[j][i];
    }
  }

  r_ECI = std::vector<double>(3, 0.0);
  v_ECI = std::vector<double>(3, 0.0);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      r_ECI[i] += R_transpose[i][j] * r_perifocal[j];
      v_ECI[i] += R_transpose[i][j] * v_perifocal[j];
    }
  }

  // Print the resulting vectors r_ECI and v_ECI
  // format data output
  std::cout << "\t\t----------------------------------" << std::endl;
  std::cout << "\t\t...PRINTING INITIAL POSITION...\n" << std::endl;

  std::cout << "\t\t -- Perifocal RF --\n" << std::endl;
  std::cout << "\t\tr0_p = [" << r_perifocal[0] << ", " << r_perifocal[1]
            << ", " << r_perifocal[2] << "]" << std::endl;
  std::cout << "\t\tv0_p = [" << v_perifocal[0] << ", " << v_perifocal[1]
            << ", " << v_perifocal[2] << "]\n"
            << std::endl;

  std::cout << "\t\t -- ECI RF --\n" << std::endl;
  std::cout << "\t\tr0_ECI = [" << r_ECI[0] << ", " << r_ECI[1] << ", "
            << r_ECI[2] << "]" << std::endl;

  std::cout << "\t\tv0_ECI = [" << v_ECI[0] << ", " << v_ECI[1] << ", "
            << v_ECI[2] << "]\n"
            << std::endl;

  // Initialise
  double temp = 0.0;
  DU = 0.0;
  TU = 0.0;

  for (const auto &x : r_ECI) {
    temp += x * x;
  }
  DU = std::sqrt(temp);
  TU = std::sqrt((DU * DU * DU) / Earth_mu);

  // Initialise
  r_0 = std::vector<double>(3, 0.0);
  v_0 = std::vector<double>(3, 0.0);

  // Undimensionalise coordinates
  for (int i = 0; i < 3; i++) {
    r_0[i] = r_ECI[i] / DU;
    v_0[i] = v_ECI[i] / (DU / TU);
  }

  std::cout << "\t\t -- UNDIMENSIONALISE --\n" << std::endl;
  std::cout << "\t\tr0 = [" << r_0[0] << ", " << r_0[1] << ", " << r_0[2] << "]"
            << std::endl;

  std::cout << "\t\tv0 = [" << v_0[0] << ", " << v_0[1] << ", " << v_0[2]
            << "]\n"
            << std::endl;
  std::cout << "\t\t----------------------------------\n" << std::endl;

  x = r_0;
  x.insert(x.end(), v_0.begin(), v_0.end());
};

/**
 * @brief
 *
 */
std::vector<double> Orbit::EoM(std::vector<double> &x) {

  double R_E = Earth_Radius / SMA;
  double mu = 1.0;

  double temp1 = 0.0;
  for (int i = 0; i < 3; ++i) {
    temp1 += x[i] * x[i];
  }
  double r_norm = sqrt(temp1);

  // calculate and aggregate the acceleration from gravity and J2
  a_total = std::vector<double>(3, 0.0);
  a_gravity = std::vector<double>(3, 0.0);
  a_oblateness = std::vector<double>(3, 0.0);

  double common_term = -(3 * J2 * mu * R_E * R_E) / (2 * std::pow(r_norm, 5));
  double z_squared_term = (5 * std::pow(x[2], 2)) / std::pow(r_norm, 2);

  a_oblateness[0] = common_term * x[0] * (1 - z_squared_term);
  a_oblateness[1] = common_term * x[1] * (1 - z_squared_term);
  a_oblateness[2] = common_term * x[2] * (3 - z_squared_term);

  for (int i = 0; i < 3; ++i) {
    a_gravity[i] = -mu * x[i] / pow(r_norm, 3);
    a_total[i] = a_gravity[i] + a_oblateness[i];
  }

  x_dot = std::vector<double>(6, 0.0);
  for (int i = 0; i < 3; ++i) {
    x_dot[i] = x[i + 3];
    x_dot[i + 3] = a_total[i];
  }

  return x_dot;
};

/**
 * @brief
 *
 *
 */
void Orbit::RungeKutta45(double dt, double T, std::vector<double> &x) {

  std::vector<double> k1 = std::vector<double>(6, 0.0);
  std::vector<double> k2 = std::vector<double>(6, 0.0);
  std::vector<double> k3 = std::vector<double>(6, 0.0);
  std::vector<double> k4 = std::vector<double>(6, 0.0);
  std::vector<double> temp = std::vector<double>(6, 0.0);

  double timestep = T / dt;

  // Display Information on output text file
  std::cout << "\t\tWriting output to file 'output.txt'." << std::endl;
  std::ofstream vOut("output.txt", std::ios::out | std::ios::trunc);

  for (int t = 0; t < timestep; t++) {

    k1 = Orbit::EoM(x);

    for (int i = 0; i < 6; i++)
      temp[i] = x[i] + (dt * k1[i]) / 2;

    k2 = Orbit::EoM(temp);

    for (int i = 0; i < 6; i++)
      temp[i] = x[i] + (dt * k2[i]) / 2;

    k3 = Orbit::EoM(temp);

    for (int i = 0; i < 6; i++)
      temp[i] = x[i] + dt * k3[i];

    k4 = Orbit::EoM(temp);

    for (int i = 0; i < 6; i++) {
      x[i] += (1.0 / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * dt;
    }

    if (t % 628 == 0) {
      std::cout << "\t\tTime-step: " << std::setw(6) << t << "/" << timestep
                << " (" << std::setw(2) << 100 * t / timestep << "%)"
                << std::endl;
    }

    // Checking that file opened successfully.
    if (vOut.is_open()) {
      vOut << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4]
           << " " << x[5] << std::endl;
    } else {
      std::cout << "Did not open vOut successfully!" << std::endl;
    }
  }

  // Finish writing and close file
  vOut.close();
  std::cout << "\t\tFinished writing to file." << std::endl;
}

/**
 * @brief
 *
 */
void Orbit::integrate() {

  // define target time and dt
  T = 10 * 2 * M_PI * sqrt(pow(SMA, 3.0) / Earth_mu) / TU;
  double dt = 0.01;

  // Initialise r_final and v_final vectors
  r_final = std::vector<double>(3, 0.0);
  v_final = std::vector<double>(3, 0.0);

  // Begine RK45 Integration
  Orbit::RungeKutta45(dt, T, x);

  // Dimensionalise r and v vectors
  for (int i = 0; i < 3; ++i) {
    r_final[i] = x[i] * DU;
    v_final[i] = x[i + 3] * (DU / TU);
  }

  std::cout << "\n\t\trf = [" << r_final[0] << ", " << r_final[1] << ", "
            << r_final[2] << "]" << std::endl;

  std::cout << "\t\tvf = [" << v_final[0] << ", " << v_final[1] << ", "
            << v_final[2] << "]\n"
            << std::endl;
  std::cout << "\t\t----------------------------------" << std::endl;
}

/**
 * @brief Perform 3x3 matrix multiplication
 *
 * @return std::vector<std::vector<double> >
 */
std::vector<std::vector<double>>
Orbit::matrixmultiply(std::vector<std::vector<double>> A,
                      std::vector<std::vector<double>> B) {
  std::vector<std::vector<double>> result(A.size(),
                                          std::vector<double>(B[0].size(), 0));
  if (A[0].size() == B.size()) {
    for (int i = 0; i < A.size(); i++) {
      for (int j = 0; j < B[0].size(); j++) {
        for (int k = 0; k < A[0].size(); k++)
          result[i][j] += A[i][k] * B[k][j];
      }
    }
  }

  return result;
};
