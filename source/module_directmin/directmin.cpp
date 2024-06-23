#include "directmin.h"
#include <iostream>

namespace ModuleDirectMin{


    DirectMin::DirectMin(Input & inp, UnitCell & cell)
    {
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