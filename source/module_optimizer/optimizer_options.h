#ifndef OPTIMIZER_OPTIONS_H
#define OPTIMIZER_OPTIONS_H


#include <string>

// default using armadillo for linear algebra operations
// in future implement our own or other library
// #ifndef DIRECTMIN_USE_ARMA
// #define DIRECTMIN_USE_ARMA
// #endif

namespace Module_Optimizer
{
    class Options
    {
        public:

        std::string choice; // Optimizer choice, either trust-region (tr) or line-search (ls) for now, tr is for later though.
        std::string retraction_type; // Retraction for the Stiefel manifold
        std::string vectransport_type; // Vector transport for the Stiefel manifold

        // 'test': test problems
        // 'ks': Kohn-Sham energy functional 
        // 'ks-sic' : self-interaction-corrected Kohn-Sham energy functional
        // 'hybrid': hybrid Kohn-Sham energy functional
        // 'rdmft': reduced density matrix energy functional
        std::string obj_type; 

        Options();
        virtual void print_info() ;
    };
};



#endif // OPTIMIZER_OPTIONS_H