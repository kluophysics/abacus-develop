#include "linesearch_helper.h"

#include <algorithm>    // std::swap
#include <limits>       // std::numeric_limits

namespace ModuleDirectMin
{

    double find_quadratic_minimizer(double alpha_lo, double phi_lo, double dphi_lo, double alpha_up, double phi_up) 
    {
        return alpha_lo + (((alpha_up - alpha_lo)*(alpha_up - alpha_lo))*dphi_lo)/(2.*(phi_lo - phi_up + (alpha_up - alpha_lo)*dphi_lo));
    }

    double find_cubic_minimizer(double alpha_lo, double phi_lo, double dphi_lo, double alpha_up, double phi_up, double dphi_up)
    {        
        if (alpha_lo>alpha_up) 
        {
            std::swap(alpha_lo, alpha_up);
            std::swap(phi_lo,   phi_up);
            std::swap(dphi_lo, dphi_up);
        }
        double d1 = 3.*(phi_lo - phi_up)/(alpha_up - alpha_lo) + dphi_lo + dphi_up;
        double D = d1*d1 - dphi_lo*dphi_up;
        if (D<=0) return std::numeric_limits<double>::max(); // no minumum in the interval, +inf here because of the linesearch nature
        double d2 = std::sqrt(D); // this code assumes alpha_lo<alpha_up; negate this value if alpha_up<alpha_lo.
        return alpha_up - ((alpha_up - alpha_lo)*(dphi_up + d2 - d1))/(dphi_up - dphi_lo + 2.*d2);
    }


    double find_quadratic_minimizer(double alpha_lo, double dphi_lo, double alpha_up, double dphi_up) 
    {
        return alpha_up + ((alpha_up - alpha_lo)*dphi_up)/(dphi_lo - dphi_up);
    };

    bool sufficient_decrease(double step_size, double dphi_lo, double phi_lo, double phi_up, double c1) 
    {
        return phi_up  <= phi_lo + c1 * step_size * dphi_lo;
    }

    bool curvature_condition(double dphi_lo, double dphi_up, double c2) {
        return std::abs(dphi_up) <= c2*std::abs(dphi_lo);
    }

}