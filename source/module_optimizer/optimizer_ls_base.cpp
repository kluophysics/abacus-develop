#include "optimizer_ls_base.h"

#include <iomanip>
#include <list>


#include "manifolds/manifold.h"
#include "module_base/timer.h"
#include "module_base/constants.h"

namespace Module_Optimizer
{
    // using ManifoldPoint = Manifold::ManifoldPoint;
    // using ManifoldVector = Manifold::ManifoldVector;

	void OptimizerLSBase::optimize()
	{
		pre_funs.clear();

		// std::cout << "--------------------------------"
		// << "inside " << "OptimizerLSBase::optimize()" << std::endl;
		ManifoldPoint xTemp;
        ManifoldVector  gfTemp;
		// x1.brief_print();
		// f1.brief_print();
		f1 = prob->objective_function(x1); nf ++; // f value at x1
		f2 = f1;
		gf1 = prob->rie_grad(x1); ng ++;// grad value at x1
		// gf1 = prob->grad(x1); ng ++;// grad value at x1

		// std::cout << "norm of x1" << x1.norm() << std::endl;

		// gf1.brief_print("gf1:");
		// ngf0= sqrt(metric(x1, gf1, gf1)); // norm of the grad 
		ngf0= sqrt(mani->metric(x1, gf1, gf1)); // norm of the grad 

		ngf1 = ngf0;  ngf2 = ngf1;
		new_slope = 0.0;


		step_size_old = 1.0;
		initial_slope_pre = 0.0;

		iter = 0;
		// std::cout << "max_iterations = " << max_iterations << std::endl;
		while (iter < max_iterations) 
		{
			// std::cout << "--------------------------------"
			// std::cout << "iter/maxiter: " << iter << "/" << max_iterations << std::endl;
            ModuleBase::timer::tick("OptimizerLSBase", "iteration");

			get_search_direction(); // obtain search direction d1;

			initial_slope = mani->metric(x1, gf1, d1);


			init_step_guess();

			// if(iter == 1)
			// 	step_size = initial_step_size;
			// else
			// {
			// step_size = 1.01*2.0*(f1-f_previous) / initial_slope; //
			// 	step_size = (step_size > 1) ? 1 : step_size;
			// 	step_size = (step_size < std::numeric_limits<double>::epsilon()) ? 1 : step_size;
			// }	

			initial_length = step_size;

			// //*If accurate enough, then a fixed step_size is chosen.*//
			// if (ngf1 /  (ngf0 + tolerance) < accuracy)
			// {
			// 	step_size = (final_step_size > 0) ? final_step_size : initial_length;
			// 	f2 = phi();
			// 	gf2 = prob->rie_grad(x2); ng++;
			// }
			// else
			// {
			// 	do_line_search();
			// }


			do_line_search();
			ngf2 = sqrt(mani->metric(x2, gf2, gf2));

			// step_size = 0.0001;
			update_data();
			// x2.brief_print("x2:");
			// x1.brief_print("x1:");
			// iter ++;
			// update_data();
			// std::cout << "norm of gf2" << ngf2 << std::endl;
			// if(ngf2 < 0 )
			// {
			// 	continue; // skip rest.
			// }
			if (std::isnan(f2) || std::isinf(f2)) /*Stop when got a nan or inf*/
			{
				printf("New function value is either nan or inf. Stop!\n");
				break;
			}

			if(0)
			{
				std::cout << "ngf1 " <<  ngf1
						  << " ngf2 " << ngf2 
						  << " f1 " << f1 
						  << " f2 " << f2 
						  << " gtol " << gtol
						  << " ftol " << ftol 
						  << std::endl;
			}
			


			// if( (abs(ngf2/ngf1 - 1) < gtol) || (abs(f1 - f2) < ftol ) )
			// {
			// 	// std::cout << (fabs(ngf1 - ngf2) < gtol) << " " << (fabs(f1 - f2) < ftol) ;
			// 	std::cout << "breaking out while loop " << std::endl;
			// 	break;
			// }
			print_info();

			// if( (fabs(ngf2/ngf1 - 1) < gtol ) || (fabs(f1 - pre_funs.back()) < ftol  ) )
			// if( iter > 1 && ( (fabs(ngf2/ngf1 - 1) < gtol ) || (fabs(f1 -f2) < ftol  ))  )
			// if( iter > 1 && ( (fabs(ngf2/ngf1 - 1) < gtol ) || (f1 - f2) < ftol  )  )
			// if( iter > 1 &&  (fabs(ngf1) < gtol ) || (fabs(f1 -f2) < ftol  ) )
			if( iter > 0 &&  (fabs(ngf1) < gtol ) || f1 - f2 < ftol   )
			// if( iter > 1 &&  (fabs(ngf1) < gtol )  )
			// if( iter > 1 &&  (fabs(ngf1) < gtol ) && (fabs(f1 -f2) < ftol  ) )
			// if( (fabs(ngf1) < gtol ) || false  )
			{
				if (1)
				{
					if (fabs(ngf1) < gtol )
						std::cout << "gtol fulfilled: " << "gtol = " << gtol  << std::endl;
					else if (f1 -f2 < ftol  )
						std::cout << "ftol fulfilled: " << "ftol = " << ftol << std::endl;
				}

				// printf("ngf2/ngf1 = %e, gtol = %e, ftol = %e, \n", fabs(ngf2/ngf1 - 1), gtol, ftol );
				// std::cout << (fabs(ngf1 - ngf2) < gtol) << " " << (fabs(f1 - f2) < ftol) ;
				std::cout << "breaking out while loop " << std::endl;
				break;
			}

			// printf("%03d/%03d %10.6f %10.4e %10.4e %-10.4e %10.4e %10.4e\n", iter, max_iterations, step_size,
			// 	ngf1, fabs(ngf2-ngf1), f2 - f1, arma::norm(x1.t()*x1)-1, f1
			// 	);
			xTemp = x1; x1 = x2; x2 = xTemp;
			gfTemp = gf1; gf1 = gf2; gf2 = gfTemp;

			pre_funs.push_front(f1);

			// f_previous = f1;
			f1 = f2;
			ngf1 = ngf2;
            initial_slope_pre = initial_slope;
			// step_size_old = (step_size_old > step_size) ? step_size_old : step_size ; // store the old step size, try next time directly with this one
			step_size_old = step_size ; // store the old step size, try next time directly with this one
			iter ++;

			ModuleBase::timer::tick("OptimizerLSBase", "iteration");

		}
		OptimizerLSBase::print_info(); // summary of nf, ng, nV, nR, nH.
		
	}

	void OptimizerLSBase::print_info()
	{
		// adapt using ModuleBase later...
		// printf("%03d/%03d %6.3E %10.4E %10.4E %-10.4E %5.3f %-15.10f\n", iter, max_iterations, step_size,
		// ngf1, fabs( ngf2/ngf1 - 1), (f2 - f1)* ModuleBase::Ry_to_eV, x1.norm(), f1* ModuleBase::Ry_to_eV
		// );

		printf("%03d/%03d %6.3E %10.4E %10.4E %-10.4E  %-15.10f\n", iter, max_iterations, step_size,
		ngf1, fabs( ngf2/ngf1 - 1), (f2 - f1)* ModuleBase::Ry_to_eV, f1* ModuleBase::Ry_to_eV
		);
	}

    void OptimizerLSBase::set_default_params()
    {
        OptimizerBase::set_default_params();
		condition_type = STRONG_WOLFE;
		LS_alpha = 1e-4;
		LS_beta = 0.999;
		LS_c1 = 1e-4;
		LS_c2 = 0.9;
		initial_step_size = 1.0;
		max_step_size = 1e10;
		min_step_size = std::numeric_limits<double>::epsilon();
		final_step_size = -1.0;
		gtol = 1e-6;
		ftol = 1e-8;
		tolerance = 1e-5;
		accuracy = 1e-8;
		step_size  = 1.0;
        return ;
    }

    void OptimizerLSBase::update_params(LSOptions *opt_in)
    {
        OptimizerBase::update_params(opt_in);
		if(opt_in ->ls_condition == "swolfe" )
		{
			condition_type = STRONG_WOLFE;
		}
		else if(opt_in ->ls_condition == "armijo")
		{
			condition_type = ARMIJO;
		}
		else if(opt_in ->ls_condition == "wolfe")
		{
			condition_type = WOLFE;
		}
		else
		{
			ModuleBase::WARNING("OptimizerLSBase::update_params", "Unknown line search condition. Use Strong Wolfe condition.");
			condition_type = STRONG_WOLFE;
		}

		LS_alpha = opt_in -> ls_alpha;
		LS_beta  = opt_in -> ls_beta;
		LS_c1    = opt_in -> ls_c1;
		LS_c2    = opt_in -> ls_c2;

		initial_step_size =  opt_in -> ls_initstepsize;
		max_step_size     =  opt_in -> ls_maxstepsize;
		min_step_size  	  =  opt_in -> ls_minstepsize ;
		final_step_size   =  opt_in -> ls_finalstepsize;

		gtol = opt_in -> ls_gtol;
		ftol = opt_in -> ls_ftol;


        return ;
    }

    void OptimizerLSBase::do_line_search()
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
		else
		{
			ModuleBase::WARNING("OptimizerLSBase::do_line_search", "Unknown line search condition. Use Strong Wolfe condition.");
		}
        return ;
    }

    void OptimizerLSBase::init_step_guess()
    {
		double alpha = 1.0;

		if (iter == 0)
		{
			step_size = initial_step_size;
		}
		else {

	
		// alpha = 2.0*(f2 - f1)/ slope;
		alpha = (1.01*2.0) *(f1 - pre_funs.front() ) / initial_slope_pre;
		// alpha = (alpha > max_step_size) ? max_step_size : alpha;

		// step_size = (1.2*step_size_old < alpha) ? alpha : 1.2*step_size_old;
		// step_size = (step_size_old < alpha) ? alpha : step_size_old;
		// step_size = (step_size > max_step_size) ? step_size : max_step_size;
		// step_size = (alpha > 1.0 ) ? 1.0 : alpha;
		// step_size = (alpha > 1.0 ) ? 1.0 : alpha;

		step_size = alpha;
		// step_size = step_size;
		}
    }

    double OptimizerLSBase::phi()
    {
		d2 = step_size * d1;
		x2 = mani->retraction(x1, d2); nR++;
		nf ++;
		return prob->objective_function(x2); // return f(x2);
    }

    double OptimizerLSBase::dphi()
    {
		// evaluate in the Stiefel manifold, not in Eucliean space
		d2 = step_size * d1;
		// x2 = retraction(x1, d2);
		gf2 = prob->rie_grad(x2); 	ng++;
		// gf2 = prob->rie_grad(x2); ng ++;
		// ManifoldPoint diff_d2 = diff_retraction(x2, d2); nV++;

		ManifoldVector diff_d2 = mani->diff_retraction(x2, gf2, x2, gf2) ; nV++;
		return mani->metric(x2, gf2, diff_d2);
    }




	void OptimizerLSBase::zoom(
		    double x1,
        double fx1,
        double slopex1, 
        double x2, 
        double fx2)
    {
		double xdiff, xincr, xlo = x1, xhi = x2, fxlo = fx1, fxhi = fx2, xlo_slope = slopex1;
		int times = 0;
		while (1)
		{
			// can use cubic interpolation to shorten the evaluation!!!
			xdiff = (xhi - xlo);
			xincr = -xlo_slope * xdiff * xdiff / 2 / (fxhi - (fxlo + xlo_slope * xdiff));
            xincr = (xincr < xdiff * LS_c1) ? xdiff * LS_c1 : xincr;
            xincr = (xincr < xdiff * LS_c2) ? xincr : xdiff * LS_c2;
			step_size = xlo + xincr;
			times++;
			if (times >= 10)
			{
				// LSstatus = LSSM_LSERROR;
				return;
			}
			f2 = phi();
			if (f2 > f1 + LS_alpha * step_size * initial_slope || f2 >= fxlo)
			{
				xhi = step_size;
				fxhi = f2;
			}
			else
			{
				new_slope = dphi();
				if (fabs(new_slope) <= -LS_beta * initial_slope)
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
			if (step_size <= min_step_size)
			{
				// step_size = min_step_size;
				// LSstatus = LSSM_MINSTEPSIZE;
				return;
			}
		};
	}

	void OptimizerLSBase::StrongWolfe()
	{
		double previous_step_size = 0, 
		f_previous = f1, 
		newslope_previous = initial_slope;
        double tstep_size = 0;
		// LSstatus = LSSM_SUCCESS;
		while (1)
		{
			// if(verbose)
			// {
			// 	std::cout << " Inside OptimizerLSBase::StrongWolfe()" << std::endl;
			// }
			f2 = phi();  
			if (f2 > f1 + LS_alpha * step_size * initial_slope || f2 >= f_previous)
			{
				zoom(previous_step_size, f_previous, newslope_previous, step_size, f2);
				return;
			}
			new_slope = dphi(); 
			if (fabs(new_slope) <= -LS_beta * initial_slope)
			{
				return;
			}
			if (new_slope >= 0)
			{
				zoom(step_size, f2, new_slope, previous_step_size, f_previous);
				return;
			}
			previous_step_size = step_size;
			f_previous = f2;
			newslope_previous = new_slope;
			if (step_size >= max_step_size)
			{
				// LSstatus = LSSM_MAXSTEPSIZE;
				return;
			}
            tstep_size = - initial_slope * step_size * step_size / 2 / (f2 - f1 - initial_slope * step_size);
            step_size = (1.1 * step_size < tstep_size) ? tstep_size : 1.1 * step_size;
			step_size = (step_size < max_step_size) ? step_size : max_step_size;
		}
	}

	void OptimizerLSBase::Wolfe()
	{

	}


	void OptimizerLSBase::Armijo()
	{
	}
}