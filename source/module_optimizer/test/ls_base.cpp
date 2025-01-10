#include "optimizer_ls_base.h"
#include <cmath>

namespace Module_Optimizer {

    void Optimizer_LS_Base::Armijo() {
        double alpha = 1.0; // Initial step size
        double beta = 0.5;  // Step size reduction factor
        double sigma = 1e-4; // Armijo condition parameter

        ManifoldPoint x = current_point; // Current point in the manifold
        ManifoldVector grad = current_gradient; // Current gradient
        ManifoldVector direction = search_direction; // Search direction

        double f_x = objective_function(x); // Objective function value at current point
        double f_x_new;
        ManifoldPoint x_new;

        while (true) {
            x_new = manifold.retraction(x, alpha * direction); // Compute new point
            f_x_new = objective_function(x_new); // Objective function value at new point

            // Check Armijo condition
            if (f_x_new <= f_x + sigma * alpha * manifold.metric(x, grad, direction)) {
                break; // Armijo condition satisfied
            }

            alpha *= beta; // Reduce step size
        }

        // Update current point and step size
        current_point = x_new;
        step_size = alpha;
    }

    void Optimizer_LS_Base::Wolfe() {
        double alpha = 1.0; // Initial step size
        double beta = 0.5;  // Step size reduction factor
        double sigma = 1e-4; // Armijo condition parameter
        double c2 = 0.9; // Curvature condition parameter

        ManifoldPoint x = current_point; // Current point in the manifold
        ManifoldVector grad = current_gradient; // Current gradient
        ManifoldVector direction = search_direction; // Search direction

        double f_x = objective_function(x); // Objective function value at current point
        double f_x_new;
        ManifoldPoint x_new;
        ManifoldVector grad_new;

        while (true) {
            x_new = manifold.retraction(x, alpha * direction); // Compute new point
            f_x_new = objective_function(x_new); // Objective function value at new point
            grad_new = gradient(x_new); // Gradient at new point

            // Check Armijo condition
            if (f_x_new > f_x + sigma * alpha * manifold.metric(x, grad, direction)) {
                alpha *= beta; // Reduce step size
                continue;
            }

            // Check curvature condition
            if (manifold.metric(x_new, grad_new, direction) < c2 * manifold.metric(x, grad, direction)) {
                alpha *= beta; // Reduce step size
                continue;
            }

            break; // Both conditions satisfied
        }

        // Update current point and step size
        current_point = x_new;
        current_gradient = grad_new;
        step_size = alpha;
    }

    void Optimizer_LS_Base::StrongWolfe() {
        double alpha = 1.0; // Initial step size
        double beta = 0.5;  // Step size reduction factor
        double sigma = 1e-4; // Armijo condition parameter
        double c2 = 0.9; // Strong Wolfe curvature condition parameter

        ManifoldPoint x = current_point; // Current point in the manifold
        ManifoldVector grad = current_gradient; // Current gradient
        ManifoldVector direction = search_direction; // Search direction

        double f_x = objective_function(x); // Objective function value at current point
        double f_x_new;
        ManifoldPoint x_new;
        ManifoldVector grad_new;

        while (true) {
            x_new = manifold.retraction(x, alpha * direction); // Compute new point
            f_x_new = objective_function(x_new); // Objective function value at new point
            grad_new = gradient(x_new); // Gradient at new point

            // Check Armijo condition
            if (f_x_new > f_x + sigma * alpha * manifold.metric(x, grad, direction)) {
                alpha *= beta; // Reduce step size
                continue;
            }

            // Check strong Wolfe curvature condition
            if (std::abs(manifold.metric(x_new, grad_new, direction)) > c2 * std::abs(manifold.metric(x, grad, direction))) {
                alpha *= beta; // Reduce step size
                continue;
            }

            break; // Both conditions satisfied
        }

        // Update current point and step size
        current_point = x_new;
        current_gradient = grad_new;
        step_size = alpha;
    }

}