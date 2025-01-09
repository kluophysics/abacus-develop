#include "optimizer_options.h"

#include <iostream>

Module_Optimizer::Options::Options()
{
        choice = "lr";
        retraction_type = "qr";
        vectransport_type = "qr";
        // maxiter = 100;
        obj_type = "test";
}

void Module_Optimizer::Options::print_info()
{
        if(choice == "lr")
        {
            std::cout << "Using Line Search Algorithm" << std::endl; 
            std::cout << "Stiefel Manifold retraction type  is "
                << retraction_type << std::endl;
            std::cout << "vector transport type is "
                << vectransport_type << std::endl;
        }
        else if (choice == "tr")
        {
            std::cout << "Using Trust Region Algorithm" << std::endl; 
            std::cout << "Stiefel Manifold retraction type  is "
                << retraction_type << std::endl;
            std::cout << "vector transport type is "
                << vectransport_type << std::endl;
        }
        else
        {
            std::cout << "Unknown choice, only lr or tr is allowed" << std::endl;
        }
}
