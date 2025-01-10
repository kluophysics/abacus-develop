#ifndef OPTIMIZER_LS_OPTIONS_H
#define OPTIMIZER_LS_OPTIONS_H


#include <string>

#include "optimizer_options.h"


namespace Module_Optimizer
{
  
    class LSOptions : public Options
    {
    public:
        // LSOptions(Input & input);  // constructor
        ~LSOptions(); // destructor

    public:

        std::string ls_method; // DirectMin line search method
        std::string ls_condition; // DirectMin line search condition
        std::string ls_cg_algo; // DirectMin line search method cg algorithm
        std::string ls_initstep_type; // initial step type for DirectMin line search method, such as 
        double ls_alpha; // the coefficient of the Wolfe first condition
        double ls_beta;  // the coefficient of the Wolfe second condition
        double ls_c1; // the coefficient of the Armijo-Goldstein condition
        double ls_c2; // the coefficient of the Armijo-Goldstein condition

        double ls_minstepsize; //the minimum stepsize allowed in the linesearch algorithm
        double ls_maxstepsize; //the maximum stepsize allowed in the linesearch algorithm
        double ls_initstepsize; // initial stepsize at the first iteration
        double ls_finalstepsize; // final stepsize 

        double ls_gtol; // line search gradient tolerance
        double ls_ftol; // line search function value tolerance

        void print_info();
    };

};



#endif // OPTIMIZER_LS_OPTIONS_H