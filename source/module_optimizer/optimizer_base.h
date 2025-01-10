#ifndef OPTIMIZER_BASE_H
#define OPTIMIZER_BASE_H

#include <string>

#include "problems/problem.h"

#include "module_parameter/parameter.h"
#include "module_parameter/input_parameter.h"

#include "module_cell/unitcell.h"

#include "optimizer_options.h"
#include "optimizer_misc.h"
#include "optimizer_logger.h"

#include "manifolds/manifold.h"

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

    /*  Objective function type
    TEST: test problems
    TEST_EIG: test problems with eigenvalues
    KS: Kohn-Sham energy functional
    KS_SIC: self-interaction-corrected Kohn-Sham energy functional
    KS_HYBRID: hybrid Kohn-Sham energy functional
    RDMFT: reduced density matrix energy functional
    */
    enum ObjectiveType { TEST, TEST_EIG, KS, KS_SIC, KS_HYBRID, RDMFT, ObjectiveTypeLength};

    /*  Condition type for line search
    ARMIJO: The Armijo-Goldstein condition.[DS83 Algorithm A6.3.1] combined with nonmontone line search
    WOLFE: The weak Wolfe condition [DS83 Algorithm A6.3.1mod]
    STRONG_WOLFE: The strong Wolfe condition [NW06 Algorithm 3.5]
    */
    enum ConditionType {ARMIJO, WOLFE,  STRONG_WOLFE, ConditionTypeLength};
                    


class OptimizerBase
{
public: 
    virtual ~OptimizerBase() = 0;

    // initialize( the problem involved
    virtual void initialize() =0 ;

    // this is to interface with the electronic structure calculation
    virtual void initialize(const Input_para & inp, UnitCell & cell) =0 ;

    // virtual void initialize(Problem *prob_in, Logger * logger_in, Options * options_in);

    // set default parameters
    virtual void set_default_params() = 0;
    // update parameters with Options pointer options_in
    virtual void update_params(Options* options_in)=0;
    // the working function for the optimization
    virtual void optimize() = 0;

protected:
    std::string name; // name of the optimizer
    
    std::size_t max_iterations; // max number of iterations
    std::size_t min_iterations; // max number of iterations

    StopCriterion stop_criterion; // stopping criterion
    VerbosityLevel verbosity; // verbosity level
    OptimizerType optimizer_type;    // optimization type, either LS or TR


    // Problem<std::complex<double>> * prob_cplx;
    // Problem<double> * prob;
    Problem * prob;
    Logger * logger;
    Options * options;

    Manifold * mani;

    // ManifoldPoint x1, x2;

    
};

// void optimize();


}


#endif  // OPTIMIZER_BASE_H