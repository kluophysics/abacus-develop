#include "directmin.h"
#include <iostream>

#include "module_base/timer.h"

namespace ModuleDirectMin{


    DirectMin::DirectMin(Input & inp, UnitCell & cell)
    {

        // ModuleBase::timer::tick("ModuleDirectMin", "DirectMin::DirectMin");

        // initialize LineSearchOptions
        if(inp.directmin_choice == "ls")
        {
            ls_opts = new LineSearchOptions(inp);
        }
        else if (inp.directmin_choice == "tr")
        {
            std::cout << "Warning: trust region has not been implemented" << std::endl; 
            // exit();
        }
        else
        {
            std::cout << "ERROR: Unknown dm_choice " << inp.directmin_choice << std::endl;
        }


        // next determine which line-search algorithm to use
        // currently works only for the Riemannian conjugate gradient (RCG) method
        // and the Riemannian steepest descent (RSD)

        if(ls_opts->ls_method == "sd")
        {
            // optimizer = new RCG(prob, ls_opts);
        }
        else if(ls_opts->ls_method == "cg")
        {
            // optimizer = new RCG(prob, ls_opts);
        }


        // choose objective type 
        if (ls_opts->obj_type == "test")
        {
            // work for test problem
        }
        else if(ls_opts-> obj_type == "ks")
        {
            // work for Kohn-Sham problem 

            // std::cout << "Warning: KS not implemented" << std::endl;
            //prob = new Problem_KS(inp, cell);
        }
        else if (ls_opts -> obj_type == "rdmft")
        {
            // prob = NULL;
            std::cout << "Warning: RDMFT not implemented" << std::endl;
        }
        else 
        {   
            // prob = NULL;
            std::cout << "ERROR: unknown objective functional type" << std::endl;
        }


        // check which 

    }
    // {
    //     // std::cout << "dm_choice: " << inp.dm_choice << std::endl;

    //     // input = inp;
    //     if(inp.dm_choice == "ls")
    //     {
    //         ls_opts = new LineSearchOptions(inp);
    //     }
    //     else if (inp.dm_choice == "tr")
    //     {
    //         std::cout << "Warning: trust region has not been implemented" << std::endl; 
    //     }
    //     else
    //     {
    //         std::cout << "ERROR: Unknown dm_choice " << inp.dm_choice << std::endl;
    //     }
    // }

    void DirectMin::initialize(Input & inp, UnitCell & cell)
    {

    //     // define the problem first!!!
    //     if(ls_opts->obj_type == "test")
    //     {
    //         prob = new Problem_Test();
    //     }
    //     if(ls_opts->obj_type == "test_eig")
    //     {
    //         prob = new Problem_Eig();
    //     }
    //     else if(ls_opts-> obj_type == "ks")
    //     {
    //         // std::cout << "Warning: KS not implemented" << std::endl;
    //         prob = new Problem_KS(inp, cell);
    //     }
    //     else if (ls_opts -> obj_type == "rdmft")
    //     {
    //         prob = NULL;
    //         std::cout << "Warning: RDMFT not implemented" << std::endl;
    //     }
    //     else 
    //     {   
    //         prob = NULL;
    //         std::cout << "ERROR: unknown functional type" << std::endl;
    //     }

    //     // determin the algorithm used
    //     if(ls_opts->ls_method == "cg")
    //     {
    //         optimizer = new RCG(prob, ls_opts);
    //     }
    //     else if (ls_opts->ls_method == "sd")
    //     {
    //         optimizer = new RSD(prob, ls_opts);
    //     }

    //     else if (ls_opts->ls_method == "bfgs")
    //     {
    //         // std::cout << "WARNING: Not implemented yet" << std::endl;
    //         optimizer = new RBFGS(prob, ls_opts);
    //     }
    //     else
    //     {
    //         std::cout << "WARNING: Not implemented yet" << std::endl;
    //     }

    };

    void DirectMin::run( int istep, UnitCell &cell)
    {
    //     optimizer->optimize();
    };

};