#include "linesearch_base.h"

#include "linesearch_helper.h"

namespace ModuleDirectMin
{

    LineSearchBase::LineSearchBase()
    {
        method_name = "Unspecified";
        verbose = false;

    }

    void LineSearchBase::evaluate_phi_and_dphi(double *phi, double *dphi)
    {
        eta2 = eta1 * step_size;
    }


    void LineSearchBase::zoom(
        double x1, double fx1, double slope1, 
        double x2, double fx2, double slope2)
    {
		double xdiff, xincr, xlo = x1, xhi = x2, 
                fxlo = fx1, fxhi = fx2, xlo_slope = slope1, xhi_slope = slope2;
		int times = 0;
		while (1)
		{
			// can use cubic interpolation to shorten the evaluation!!!
			// xdiff = (xhi - xlo);
			// xincr = -xlo_slope * xdiff * xdiff / 2 / (fxhi - (fxlo + xlo_slope * xdiff));
            // xincr = (xincr < xdiff * ls_options->ls_c1) ? xdiff * ls_options->ls_c1 : xincr;
            // xincr = (xincr < xdiff * ls_options->ls_c2) ? xincr : xdiff * ls_options->ls_c2;
			// step_size = xlo + xincr;

            step_size = find_cubic_minimizer(xlo, fxlo, xlo_slope, xhi, fxhi, xhi_slope);

			times++;
			if (times >= 10)
			{
				// LSstatus = LSSM_LSERROR;
				return;
			}
            evaluate_phi_and_dphi(&f2, &new_slope);

			// f2 = phi();
			// if (f2 > f1 + ls_options->ls_alpha * step_size * initial_slope || f2 >= fxlo)
			if ( !sufficient_decrease(step_size, initial_slope, f1, f2, ls_options->ls_alpha)  || f2 >= fxlo )
			{
				xhi = step_size;
				fxhi = f2;
			}
			else
			{
				// new_slope = dphi();
				// if (fabs(new_slope) <= -ls_options->ls_beta * initial_slope)
				if ( curvature_condition(initial_slope, new_slope, ls_options->ls_beta) )
				{
					return;
				}
				if (new_slope * (xhi - xlo) >= 0)
				{
					xhi = xlo;
					fxhi = fxlo;
				}
				xlo = step_size;
				fxlo = f2;
				xlo_slope = new_slope;
			}
			if (step_size <= ls_options->ls_minstepsize)
			{
				// step_size = ls_options->ls_minstepsize;
				// LSstatus = LSSM_MINSTEPSIZE;
				return;
			}
		}
	}





}


