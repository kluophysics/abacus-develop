#include <armadillo>
#include <iostream>
#include <cmath>

// Function prototype for a simple objective function
double objective(const arma::vec& x);
arma::vec gradient(const arma::vec& x);

// Moré-Thuente line search implementation
bool moreThuenteLineSearch(const arma::vec& x0, const arma::vec& direction, 
                           arma::vec& x_new, double& alpha, 
                           double c1 = 1e-4, double c2 = 0.9, double alpha_max = 10.0) {
    double alpha_min = 0.0;
    alpha = 1.0; // Initial step size
    double phi0 = objective(x0);
    arma::vec grad0 = gradient(x0);
    double dphi0 = arma::dot(grad0, direction);

    if (dphi0 >= 0) {
        std::cerr << "The search direction is not a descent direction." << std::endl;
        return false;
    }

    while (true) {
        arma::vec x_temp = x0 + alpha * direction;
        double phi = objective(x_temp);
        arma::vec grad_temp = gradient(x_temp);
        double dphi = arma::dot(grad_temp, direction);

        // Wolfe conditions
        if (phi > phi0 + c1 * alpha * dphi0 || (alpha > alpha_min && phi >= objective(x0 + alpha_min * direction))) {
            alpha_max = alpha;
        } else if (std::abs(dphi) <= -c2 * dphi0) {
            x_new = x_temp;
            return true;
        } else if (dphi >= 0) {
            alpha_max = alpha;
        } else {
            alpha_min = alpha;
        }

        if (alpha_max - alpha_min < 1e-6) {
            std::cerr << "Line search did not converge." << std::endl;
            return false;
        }

        // Bisection to update alpha
        alpha = (alpha_min + alpha_max) / 2.0;
    }
}

int main() {
    // Example: quadratic function
    arma::vec x0 = {2.0, 2.0};  // Initial guess
    arma::vec direction = {-1.0, -1.0};  // Descent direction (e.g., gradient descent step)
    arma::vec x_new;
    double alpha;

    if (moreThuenteLineSearch(x0, direction, x_new, alpha)) {
        std::cout << "Line search successful!" << std::endl;
        std::cout << "Step size: " << alpha << std::endl;
        std::cout << "New position: " << x_new.t();
    } else {
        std::cerr << "Line search failed." << std::endl;
    }

    return 0;
}

// Quadratic objective function example
double objective(const arma::vec& x) {
    return arma::dot(x, x);  // f(x) = ||x||^2
}

// Gradient of the objective function
arma::vec gradient(const arma::vec& x) {
    return 2.0 * x;  // grad(f(x)) = 2x
}
