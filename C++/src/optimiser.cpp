#include <iostream>
#include <vector>
#include <nlopt.h>
#include </usr/local/include/cppad/ipopt/solve.hpp>
#include </opt/homebrew/opt/eigen/include/eigen3/Eigen/Dense>


// Dynamics model (ODE)
template <typename Scalar>
void lowThrustDynamics(Scalar t, const Scalar& x, const Scalar& u, Scalar& x_dot) {
    // Constants
    Scalar G = 6.67430e-11;  // Gravitational constant
    Scalar M = 5.972e24;     // Earth mass
    Scalar R = 6371e3;       // Earth radius

    // State derivatives
    x_dot = x[2];  // x_dot = v
    x_dot[1] = x[3];  // y_dot = w
    x_dot[2] = u[0] / x[0] - G * M / (x[0] * x[0]);  // v_dot = T/m - GM/r^2
    x_dot[3] = u[1] / x[0];  // w_dot = L / m
}

// Cost function to be minimized
template <typename Scalar>
Scalar trajectoryCost(const CppAD::vector<Scalar>& X, const CppAD::vector<Scalar>& U) {
    // Penalize deviation from the desired circular orbit velocity
    Scalar velocityCost = X[2] - 7000.0;

    // Penalize high thrust
    Scalar thrustCost = CppAD::abs(U[0]) + CppAD::abs(U[1]);

    // Combine costs
    return velocityCost + 0.01 * thrustCost;
}

// Define the IPOPT interface for optimization
template <typename Scalar>
class TrajectoryOptimization : public CppAD::ipopt::base<Scalar> {
public:
    // IPOPT implementation of the objective function
    bool eval_f(size_t, const CppAD::vector<Scalar>& X, const CppAD::vector<Scalar>& U, Scalar& f) override {
        f = trajectoryCost(X, U);
        return true;
    }

    // IPOPT implementation of the constraints
    bool eval_g(size_t, const CppAD::vector<Scalar>& X, const CppAD::vector<Scalar>& U, CppAD::vector<Scalar>& g) override {
        // Dynamics constraints
        Scalar x_dot[4];
        lowThrustDynamics(0.0, X, U, x_dot);
        for (int i = 0; i < 4; ++i) {
            g[i] = x_dot[i];
        }

        return true;
    }

    // Specify the structure of the Hessian
    bool eval_h(size_t, const Scalar&, const CppAD::vector<Scalar>&, Scalar, const CppAD::vector<Scalar>& U,
                CppAD::vector<Scalar>& Hessian) override {
        // The Hessian is often not needed for this type of problem, so it's left unimplemented here.
        return false;
    }
};

int main() {
    // Initial conditions
    CppAD::vector<double> X0(4), U0(2);
    X0[0] = 7000.0;  // Initial radius
    X0[1] = 0.0;     // Initial angle
    X0[2] = 0.0;     // Initial velocity
    X0[3] = 0.0;     // Initial angular velocity

    U0[0] = 0.1;  // Initial thrust (T/m)
    U0[1] = 0.0;  // Initial angular thrust (L/m)

    // Set up IPOPT options
    CppAD::ipopt::solve_result<Dvector> solution;
    CppAD::ipopt::solve<Dvector, TrajectoryOptimization>(
        "cppad_trajectory_optimization", X0, U0, options, solution);

    // Retrieve and print the optimized trajectory
    CppAD::vector<double> optimized_X = solution.x;
    CppAD::vector<double> optimized_U = solution.zl;

    std::cout << "Optimized trajectory:\n";
    for (size_t i = 0; i < optimized_X.size(); ++i) {
        std::cout << "X[" << i << "] = " << optimized_X[i] << "\n";
    }

    for (size_t i = 0; i < optimized_U.size(); ++i) {
        std::cout << "U[" << i << "] = " << optimized_U[i] << "\n";
    }

    return 0;
}