#include "linesearch_options.h"


namespace ModuleDirectMin
{
    LineSearchOptions::LineSearchOptions()
    {
        this->set_default_parameters();
    }

    LineSearchOptions::LineSearchOptions(Input& input)
    {
        // obj_type = input.directmin_obj_type;
        // choice = input.directmin_choice;
        // maxiter = input.directmin_maxiter;

        // // obj_type = "test"; // default for now, will be changed

        // retraction_type = input.directmin_retraction_type;
        // vectransport_type = input.directmin_vectransport_type;

        Options::update_from_input(input);

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

    void LineSearchOptions::set_default_parameters()
    {
        ls_method = "cg";
        ls_condition = "swolfe";
        if(ls_method == "cg")
        {
            ls_cg_algo = "dy";
        }
        else
        {
            ls_cg_algo="";
        }

        // ls_initstep_type = input.directmin_ls_initstep_type;

        ls_alpha = 0.1;
        ls_beta = 0.9;
        ls_c1 = 1e-4;
        ls_c2 = 0.9;

        ls_minstepsize = 1e-10;
        ls_maxstepsize = 1000;
        ls_initstepsize = 1.0;
        ls_finalstepsize = -1.0;
        ls_gtol = 1e-5;
        ls_ftol = 1e-7;
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