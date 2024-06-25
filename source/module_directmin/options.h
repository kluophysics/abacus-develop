#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>
#include "module_io/input.h"

namespace ModuleDirectMin
{

    enum RetractionMethod {
        QR,
        EXP,
        CAYLEY,
        POLAR,
        RetractionMethod_Length
    };
    // VectorTransport should  use the corresponding retraction method 

    enum VectorTranportType {
        PROJECTION,
        DIFFERENTIATED,
        VectorTranportType_Length
    };




    class Options
    {
    public:

        std::string choice; // DirectMin choice, either trust-region (tr) or line-search (ls) for now, tr is for later though.
        std::string retractionType; // Retraction for the Stiefel manifold
        std::string vectransport_type; // Vector transport for the Stiefel manifold

        std::string obj_type; // functional type
        // 'test': test problems
        // 'ks': Kohn-Sham energy functional 
        // 'ks-sic' : self-interaction-corrected Kohn-Sham energy functional
        // 'hybrid': hybrid Kohn-Sham energy functional
        // 'rdmft': reduced density matrix energy functional

        int maxiter; // max number of iterations for minimization 

        // Options(Input & input);  // constructor

        Options();
        Options(Input & input);  // constructor
        virtual void print_info();
        // virtual ~Options();
    };

} // end namespace ModuleDirectMin

#endif // !OPTIONS_H