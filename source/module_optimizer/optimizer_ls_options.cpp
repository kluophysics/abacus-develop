#include "optimizer_ls_options.h"

#include <iostream>

namespace Module_Optimizer
{
  
    LineSearchOptions::~LineSearchOptions()
    {
    }

    void LineSearchOptions::print_info()
    {
        std::cout << "Line search options: " << std::endl;
        std::cout << "Line search method: " << ls_method << std::endl;
        std::cout << "Line search condition: " << ls_condition << std::endl;
        std::cout << "Line search cg algorithm: " << ls_cg_algo << std::endl;
        std::cout << "Line search initial step type: " << ls_initstep_type << std::endl;
        std::cout << "Line search alpha: " << ls_alpha << std::endl;
        std::cout << "Line search beta: " << ls_beta << std::endl;
        std::cout << "Line search c1: " << ls_c1 << std::endl;
        std::cout << "Line search c2: " << ls_c2 << std::endl;
        std::cout << "Line search min step size: " << ls_minstepsize << std::endl;
        std::cout << "Line search max step size: " << ls_maxstepsize << std::endl;
        std::cout << "Line search initial step size: " << ls_initstepsize << std::endl;
        std::cout << "Line search final step size: " << ls_finalstepsize << std::endl;
        std::cout << "Line search gtol: " << ls_gtol << std::endl;
        std::cout << "Line search ftol: " << ls_ftol << std::endl;
        return ;
    }
};