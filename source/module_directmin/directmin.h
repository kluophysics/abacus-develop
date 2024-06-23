#ifndef DIRECTMIN_H
#define DIRECTMIN_H


#include <string>

// #include "options.h"
// #include "optimizer/options_ls.h"
#include "module_io/input.h"
#include "module_cell/unitcell.h"

// #include "optimizer/ls_rcg.h"
// #include "optimizer/ls_rsd.h"
// #include "optimizer/ls_rbfgs.h"
// #include "objective/problem.h"
// #include "objective/problem_test.h"
// #include "objective/problem_eig.h"
// #include "objective/problem_ks.h"

namespace ModuleDirectMin
{

    class DirectMin
    {
        public:


        // DirectMin();
        DirectMin(Input & inp, UnitCell & cell);
        // ~DirectMin();

        void initialize(Input & inp, UnitCell & cell);
        void run( int istep, UnitCell &cell); //  main function for the directmin procedure
        ;

        private:
        // Input input;
        // ModuleDirectMin::LineSearchOptions* ls_opts; // directmin ls_options
        // // ModuleDirectMin::LineSearchOptions* ls_opts; // directmin ls_options
        // ModuleDirectMin::OptimizerBase * optimizer;

        // ModuleDirectMin::Problem * prob;
        // RCG * test = new RCG();
    };


}



// class directmin
// {
//     public:

//     // int 
//     // 
//     int maxscf; // max number of iterations for scf 
//     std::string search_direction; // minization algorithm: 
//     // "sd": steepest decent
//     // "cg": conjugate gradient
//     // "nlcg": nonlinear conjugate gradient
//     // "blfgs": blfgs
//     double rate_occ; // rate of change in occupation
//     double rate_orb; // rate of change in orbital

//     public:

//     directmin();  // constructor
//     ~directmin(); // destructor

//     private:

//     void init(); // initialize(

//     void optimize(); // run the scf procedure
    

// };


#endif // !DIRECTMIN_H