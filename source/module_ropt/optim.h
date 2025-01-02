#ifndef OPTIM_H
#define OPTIM_H

#include <string>

#include "problem/problem.h"

#include "optim_types.h"

#include "optim_misc.h"

#include "optim_options.h"
#include "optim_logger.h"

namespace Module_Optim
{
    // #ifdef OPTIM_USE_ARMA
    //     #include <armadillo>
    //     using dvec = arma::vec; //
    // #else
    //     #include <vector>
    //     using dvec = std:vector<double>;
    // #endif

class Optimizer
{
public: 
    virtual ~Optimizer() = 0;
    virtual void optimize();

    // initialize( the problem involved
    virtual void initialize(Problem *prob_in, Logger * logger_in, Options * options_in);
    // update parameters with Options pointer options_in
    virtual void update_params(Options* options_in);
    // the working function for the optimization
    virtual void optimize();

protected:
    std::string name; // name of the optimizer
    unsigned long max_iterations; // max number of iterations
    unsigned long min_iterations; // minimum number of iterations

    StopCriterion stop_crit;
    VerbosityLevel verbosity;
    OptimizerType opt_type;


    Problem * prob;
    Logger * logger;
    Options * options;

    
};

// void optimize();


}


#endif 