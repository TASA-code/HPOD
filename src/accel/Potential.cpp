#include <cmath>
#include <iostream>
#include <fstream>
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>

#include "../Orbit/Orbit.h"
#include "Potential.h"
#include "EGM96_data.h"

using namespace Eigen;


namespace {
    double factorial(int n) {
        if (n == 0 || n == 1) {
            return 1.0;
        }

        double result = 1.0;
        for (int i = 2; i <= n; ++i) {
            result *= i;
        }

        return result;
    }
}



double Potential::Legendre_poly(int l, int m, double gamma) {
    double term = 1.0 / (std::pow(2, l) * factorial(l)) * std::pow((1.0 - gamma * gamma), m / 2);
    double term1 = std::pow(-1.0, l + m);
    double term2 = factorial(2 * l) / factorial(2 * (l - m));
    double term3 = std::pow(gamma * gamma - 1.0, l - m);

    for (int k = 1; k <= (l + m); ++k) {
        term3 *= (2 * k - 1);
    }

    double P = 0;
    P = term * term1 * term2 * term3;

    double normalised_term;
    if (m == 0) {
        normalised_term = std::sqrt((factorial(l - m) * (2 * l + 1)) / factorial(l + m)) / std::sqrt(2.0);
    } else {
        normalised_term = std::sqrt((factorial(l - m) * (2 * l + 1) * 2) / factorial(l + m));
    }

    double undimen_P_lm = normalised_term * P;
    return undimen_P_lm;
}



double Potential::dU_dr(Vector3d& r, double phi, double lambda){

    double r_norm = r.norm();
    double term, coefficient, sum;
    sum = 0.0;

    for (int l = 2; l <= l_max; ++l) {
        for (int m = 0; m <= l; ++m) {

            coefficient = EGM96_Clm[l][m] * cos(m * lambda) + EGM96_Slm[l][m] * sin(m * lambda);
            
            term = std::pow(R/r_norm,l) * (l+1) * Legendre_poly(l,m,sin(phi)) * coefficient;
            
            sum += term;
        }
    };
    sum *= -mu/pow(r_norm,2);

    return sum;
};


double Potential::dU_dphi(Vector3d& r, double phi, double lambda){

    double r_norm = r.norm();
    double term, coefficient, sum;
    sum = 0.0;

    for (int l = 2; l <= l_max; ++l) {
        for (int m = 0; m <= l; ++m) {

            coefficient = EGM96_Clm[l][m] * cos(m * lambda) + EGM96_Slm[l][m] * sin(m * lambda);
            
            term = std::pow(R/r_norm,l) * Legendre_poly(l,m+1,sin(phi)) - m*tan(phi)*Legendre_poly(l,m,sin(phi)) * coefficient;
            
            sum += term;
        }
    };
    sum *= mu/r_norm;

    return sum;
};


double Potential::dU_dlambda(Vector3d& r, double phi, double lambda){

    double r_norm = r.norm();
    double term, coefficient, sum;
    sum = 0.0;

    for (int l = 2; l <= l_max; ++l) {
        for (int m = 0; m <= l; ++m) {

            coefficient = EGM96_Slm[l][m] * cos(m * lambda) - EGM96_Clm[l][m] * sin(m * lambda);
            
            term = std::pow(R/r_norm,l) * m * Legendre_poly(l,m,sin(phi)) * coefficient;

            sum += term;
        }
    };

    sum *= -mu/r_norm;

    return sum;
};


Vector3d Potential::a(Vector3d& r, double phi, double lambda){
    
    double a_I, a_J, a_K;
    
    double dUdr, dUdphi, dUdlambda;
    
    dUdr = dU_dr(r, phi, lambda);
    dUdphi = dU_dphi(r, phi, lambda);
    dUdlambda = dU_dlambda(r, phi, lambda);

    
    double r_norm = r.norm();

    double common_term1 = (1/r_norm)*dUdr - (r[2] / (std::pow(r_norm,2.0) * sqrt(std::pow(r[0],2.0) + std::pow(r[1],2.0)))) * dUdphi;
    double common_term2 = (1 / (std::pow(r[0],2.0) + std::pow(r[1],2.0))) * dUdlambda;
    double common_term3 = (mu * r_norm) / std::pow(r_norm,3.0);

    a_I = common_term1 * r[0] - common_term2 * r[1] - common_term3;
    
    a_J = common_term1 * r[1] + common_term2 * r[0] - common_term3;

    a_K = (1.0 / r_norm) * dUdr * r[2] + (sqrt(std::pow(r[0], 2.0) + std::pow(r[1], 2.0)) / std::pow(r_norm,2.0))
            * dUdphi - common_term3;

    Vector3d a;
    a << a_I, a_J, a_K;

    return a;
}



