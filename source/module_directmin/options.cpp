#include "options.h"


namespace ModuleDirectMin
{

    Options::Options()
    {
        choice = "lr";
        retraction_type = RT_QF;
        vectransport_type = VT_PROJECTION;
        maxiter = 100;
        obj_type = KS;
    }

    Options::Options(Input & input)
    {
        update_from_input(input);
    }

    void Options::update_from_input(Input & input)
    {
        choice = input.directmin_choice;
        if( input.directmin_retraction_type == "qf")
            retraction_type = RT_QF;
        if(input.directmin_vectransport_type == "projection")
            vectransport_type = VT_PROJECTION;
        if(input.directmin_obj_type == "ks")
            obj_type = KS;
        maxiter = input.directmin_maxiter;
    }

    void Options::print_info()
    {
        if(choice == "lr")
        {
            std::cout << "Using Line Search Algorithm" << std::endl; 
            std::cout << "Stiefel Manifold retraction type  is "
                << retraction_type << std::endl;
            std::cout << "vector transport type is "
                << vectransport_type << std::endl;
            std::cout << "Maximum number of iterations is " << maxiter << std::endl;
            ;
        }
    }
} // namespace ModuleDirectMin