#include "optimizer_ls_base.h"


namespace Module_Optimizer
{
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