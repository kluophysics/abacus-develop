#include "optimizer_ls_base.h"


namespace Module_Optimizer
{
    void Optimizer_LS_Base::set_default_params()
    {
        return ;
    }

    void Optimizer_LS_Base::update_params(Options *opt_in)
    {
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