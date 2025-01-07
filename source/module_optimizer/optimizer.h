#ifndef DIRECTMIN_H
#define DIRECTMIN_H

#include <string>

#include "problems/problem.h"

#include "module_parameter/parameter.h"
#include "module_parameter/input_parameter.h"

#include "module_cell/unitcell.h"

#include "optimizer_options.h"
#include "optimizer_misc.h"
#include "optimizer_logger.h"

namespace Module_Optimizer
{

	// /*Specify what information will be output in the algorithm.
	// The value should be assigned to the member variable: "VerbosityLevel" */
    enum VerbosityLevel { LOW, MEDIUM, HIGH, COMPLETE, VerbosityLevelLength};

    // /*The algorithm is stopped when a value (specified by ther parameter) is less than the "Tolerance" (a member variable)
	// The value should be assigned to the member variable: "StopCriterion" and the applicable values are
	// FUNC_REL: |f_k - f_{k+1}| / max(|f_k|, 1)
    // FUNC_ABS: |f_k - f_{k+1}| 
	// GRAD_ABS: \|gf_k\|
	// GRAD_REL: \|gf_k\| / \|gf_0\|*/
    enum StopCriterion { FUNC_REL, FUNC_ABS, GRAD_ABS, GRAD_REL,  StopCriterionLength};

    // /*Provide two types of Optimization types
    // LS:  line search algorithm
    // TR: trust-region algorithm
    // */
    enum OptimizerType { LS, TR};


class DirectMin
{
public: 
    virtual ~DirectMin() = 0;

    // initialize( the problem involved
    virtual void initialize(const Input_para & inp, UnitCell & cell);

    // virtual void initialize(Problem *prob_in, Logger * logger_in, Options * options_in);

    // set default parameters
    virtual void set_default_params();
    // update parameters with Options pointer options_in
    virtual void update_params(Options* options_in);
    // the working function for the optimization
    virtual void optimize();

protected:
    // std::string name; // name of the optimizer
    // unsigned long max_iterations; // max number of iterations
    // unsigned long min_iterations; // minimum number of iterations

    // StopCriterion stop_crit;
    // VerbosityLevel verbosity;
    // OptimizerType opt_type;


    Problem * prob;
    Logger * logger;
    Options * options;

    
};

// void optimize();


}


#endif 