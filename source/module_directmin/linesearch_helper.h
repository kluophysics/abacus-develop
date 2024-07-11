#ifndef LINESEARCH_HELPER_H
#define LINESEARCH_HELPER_H

#include <cmath>

namespace ModuleDirectMin
{

    // minimizer of the quadratic function that interpolates f(alpha_lo), f'(alpha_lo), f(alpha_up) within the given interval
    double find_quadratic_minimizer(double alpha_lo, double phi_lo, double dphi_lo, double alpha_up, double phi_up);

    // minimizer of the cubic function that interpolates f(alpha_lo), f'(alpha_lo), f(alpha_up), f'(alpha_up) 
    // within the given interval
    double find_cubic_minimizer(double alpha_lo, double phi_lo, double dphi_lo, double alpha_up, double phi_up, double dphi_up);

    // minimizer of the quadratic function that interpolates f'(a), f'(b) within the given interval
    // N.B. the function itself is undetermined since we miss information like f(a) or f(b); however the minimizer is well-defined
    double find_quadratic_minimizer(double a, double ga, double b, double gb);

    // first strong Wolfe condition, see Eq. 3.7a, p.34 [NW06]
    bool sufficient_decrease(double step_size, double dphi_lo, double phi_lo, double phi_up, double c1) ;

    // second strong Wolfe condition, see Eq. 3.7b, p.34 [NW06]
    bool curvature_condition(double dphi_lo, double dphi_up, double c2);

}

#endif /* LINESEARCH_HELPER_H */