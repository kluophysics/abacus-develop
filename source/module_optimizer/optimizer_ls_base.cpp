#include "optimizer_ls_base.h"

#include "manifolds/manifold.h"
namespace Module_Optimizer
{
    using ManifoldPoint = Manifold::ManifoldPoint;
    using ManifoldVector = Manifold::ManifoldVector;

	void Optimizer_LS_Base::optimize()
	{
		pre_funs.clear();

		// std::cout << "--------------------------------"
		// << "inside " << "OptimizerLineSearchBase::optimize()" << std::endl;
		ManifoldPoint xTemp;
        ManifoldVector  gfTemp;
		// x1.brief_print();
		// f1.brief_print();
		f1 = prob->f(x1); nf ++; // f value at x1
		f2 = f1;
		gf1 = prob->rgrad(x1); ng ++;// grad value at x1
		// gf1 = prob->grad(x1); ng ++;// grad value at x1

		// std::cout << "norm of x1" << x1.norm() << std::endl;

		// gf1.brief_print("gf1:");
		// ngf0= sqrt(metric(x1, gf1, gf1)); // norm of the grad 
		ngf0= sqrt(metric(x1, gf1, gf1)); // norm of the grad 

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
            ModuleBase::timer::tick("OptimizerLineSearchBase", "iteration");

			get_search_direction(); // obtain search direction d1;

			initial_slope = metric(x1, gf1, d1);


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
			// 	gf2 = prob->rgrad(x2); ng++;
			// }
			// else
			// {
			// 	do_line_search();
			// }


			do_line_search();
			ngf2 = sqrt(metric(x2, gf2, gf2));

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

			ModuleBase::timer::tick("OptimizerLineSearchBase", "iteration");

		}
		OptimizerBase::print_info(); // summary of nf, ng, nV, nR, nH.
		
	}

    void Optimizer_LS_Base::set_default_params()
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

    void Optimizer_LS_Base::update_params(Options *opt_in)
    {
        OptimizerBase::update_params(opt_in);

        return ;
    }

    void Optimizer_LS_Base::do_line_search()
    {
        return ;
    }

    void Optimizer_LS_Base::init_step_guess()
    {
        return ;
    }

    double Optimizer_LS_Base::phi()
    {
        return 0.0;
    }

    double Optimizer_LS_Base::dphi()
    {
        return 0.0;
    }

    void Optimizer_LS_Base::StrongWolfe()
    {
        return ;
    }

    void Optimizer_LS_Base::Armijo()
    {
        return ;
    }

    void Optimizer_LS_Base::Wolfe()
    {
        return ;
    }

    void Optimizer_LS_Base::print_info()
    {
        return ;
    }
}