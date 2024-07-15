#include "linesearch_base.h"

#include "linesearch_helper.h"

namespace ModuleDirectMin
{

    LineSearchBase::LineSearchBase()
    {
        this->set_default_parameters();

    }
    // LineSearchBase::LineSearchBase(const Problem * prob, LineSearchOptions * ls_opt_in)
    // {
    //     this->set_default_parameters();
        
    // }

    void LineSearchBase::initialize(Problem *prob_in)
    {
        prob = prob_in;
        

		eta1 = prob->get_domain()->
        // x1 = prob->X0;
    }

    void LineSearchBase::set_default_parameters()
    {
        // parameters for linesearch conditions, either Strong Wolfe conditions,
        // or Armijo conditions
        // or Weak Wolfe conditions
        if(ls_options->ls_condition == "swolfe")
        {
            condition_type = STRONG_WOLFE;
        }
        else if(ls_options->ls_condition == "armijo")
        {
            condition_type = ARMIJO;
        }
        else if(ls_options->ls_condition == "wolfe")
        {
            condition_type = WOLFE;
        }

        // set default parameters for ls_options
        ls_options = new LineSearchOptions();

        LS_status = SUCCESS;

        method_name = "Unspecified";
        verbose = false;
    }


    void LineSearchBase::do_line_search()
	{
		if (condition_type == STRONG_WOLFE)
		{
			StrongWolfe();
		}
		else if (condition_type == WOLFE)
		{
			Wolfe();
		}
		else if (condition_type == ARMIJO)
		{
			Armijo();
		}
	}


    void LineSearchBase::evaluate_phi_and_dphi(double *phi, double *dphi)
    {
        eta2 = eta1 * step_size;
    }

    void LineSearchBase::StrongWolfe()
    {
		double step_size_previous = 0, 
                f_previous = f1, 
                slope_previous = initial_slope;
        double tstep_size = 0;
		LS_status = SUCCESS;

        // double alpha_max = 1.0;

		while (1)
		{
			if(verbose)
			{
				std::cout << " Inside OptimizerLineSearchBase::StrongWolfe()" << std::endl;
			}
            evaluate_phi_and_dphi(&f2, &new_slope);

            
            // sufficient_decrease is not satisfied, do zoom and stop
            // α∗ ←zoom(αi−1 , αi ) and stop;
            if ( !sufficient_decrease(step_size, initial_slope, f1, f2, ls_options->ls_c1)  || f2 >= f_previous )
			{
                zoom(step_size_previous, f_previous, slope_previous, step_size, f2, new_slope);
				return;
			}

            // if curvature is satisfied |φ' (αi )| ≤ −c2 φ' (0)
            //  set α∗ ← αi and stop;
            if ( curvature_condition(initial_slope, new_slope, ls_options->ls_c2) )
            {
                return;
            }


            // if φ'(αi ) ≥ 0, 
            // set α∗ ←zoom(αi , αi−1 ) and stop;
			if (new_slope >= 0)
			{
				// zoom(step_size, f2, new_slope, step_size_previous, f_previous);
                zoom(step_size_previous, f_previous, slope_previous, step_size, f2, new_slope);
				return;
			}
			step_size_previous = step_size;
			f_previous = f2;
			slope_previous = new_slope;

			if (step_size >= ls_options->ls_maxstepsize )
			{
				LS_status = MAXSTEPSIZE;
				return;
			}

            //Choose αi+1 ∈ (αi , αmax );
            // alpha_max = std::max(step_size, alpha_max);
            // step_size = 
            tstep_size = - initial_slope * step_size * step_size / 2 / (f2 - f1 - initial_slope * step_size);
            step_size = (1.1 * step_size < tstep_size) ? tstep_size : 1.1 * step_size;
			step_size = (step_size < ls_options->ls_maxstepsize) ? step_size : ls_options->ls_maxstepsize;
		}
	};


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
            step_size = find_cubic_minimizer(xlo, fxlo, xlo_slope, xhi, fxhi, xhi_slope);

			times++;
			if (times >= 10)
			{
				LS_status = LSERROR;
				return;
			}
            evaluate_phi_and_dphi(&f2, &new_slope);

			// f2 = phi();
			// if (f2 > f1 + ls_options->ls_c1 * step_size * initial_slope || f2 >= fxlo)
			if ( !sufficient_decrease(step_size, initial_slope, f1, f2, ls_options->ls_c1)  || f2 >= fxlo )
			{
				xhi = step_size;
				fxhi = f2;
			}
			else
			{
				// new_slope = dphi();
				// if (fabs(new_slope) <= -ls_options->ls_c2 * initial_slope)
				if ( curvature_condition(initial_slope, new_slope, ls_options->ls_c2) )
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
				LS_status = MINSTEPSIZE;
				return;
			}
		}
	}





}


