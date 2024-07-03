#include "options.h"


namespace ModuleDirectMin
{

    Options::Options()
    {
        choice = "lr";
        retraction_type = "qr";
        vectransport_type = "qr";
        maxiter = 100;
        obj_type = "test";
    }

    Options::Options(Input & input)
    {
        choice = input.directmin_choice;
        retraction_type = input.directmin_retraction_type;
        vectransport_type = input.directmin_vectransport_type;
        maxiter = input.directmin_maxiter;
        obj_type = input.directmin_obj_type;
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