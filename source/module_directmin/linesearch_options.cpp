#include "linesearch_options.h"


namespace ModuleDirectMin
{

    LineSearchOptions::LineSearchOptions(Input& input)
    {
        obj_type = input.directmin_obj_type;
        choice = input.directmin_choice;
        maxiter = input.directmin_maxiter;

        // obj_type = "test"; // default for now, will be changed

        retraction_type = input.directmin_retraction_type;
        vectransport_type = input.directmin_vectransport_type;

        // Options(input);
        ls_method = input.directmin_ls_method;
        ls_condition = input.directmin_ls_algo;
        ls_cg_algo = input.directmin_ls_cg_algo;

        ls_initstep_type = input.directmin_ls_initstep_type;

        ls_alpha = input.directmin_ls_alpha;
        ls_beta = input.directmin_ls_beta;
        ls_c1 = input.directmin_ls_c1;
        ls_c2 = input.directmin_ls_c2;

        ls_minstepsize = input.directmin_ls_minstepsize;
        ls_maxstepsize = input.directmin_ls_maxstepsize;
        ls_initstepsize = input.directmin_ls_initstepsize;
        ls_finalstepsize = input.directmin_ls_finalstepsize;
        ls_gtol = input.directmin_ls_gtol;
        ls_ftol = input.directmin_ls_ftol;

    }

    LineSearchOptions::~LineSearchOptions()
    {
        ;
    }

    void LineSearchOptions::print_info()
    {
        Options::print_info();

        std::cout << "  direct minimization " << ls_method  << std::endl;
        std::cout << "  choice: " << choice << std::endl;
        std::cout << "  maxiter: " << maxiter << std::endl;
        std::cout << "  functional: " << obj_type << std::endl;
        std::cout << "  retraction: " << retraction_type << std::endl;
        std::cout << "  vectransport: " << vectransport_type << std::endl;
        std::cout << "  linesearch: " << ls_method << std::endl;
        std::cout << "  linesearch condition: " << ls_condition << std::endl;
        std::cout << "  linesearch cg algorithm: " << ls_cg_algo << std::endl;
        std::cout << "  linesearch init step type: " << ls_initstep_type << std::endl;
    }
}