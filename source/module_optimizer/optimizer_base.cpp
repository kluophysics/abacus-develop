#include "optimizer_base.h"

namespace Module_Optimizer
{
    OptimizerBase::~OptimizerBase() = default;

    void OptimizerBase::set_default_params()
    {
        name = "OptimizerBase";

        max_iterations = 200;
        min_iterations = 0;

        verbosity = VerbosityLevel::MEDIUM;
        stop_criterion = StopCriterion::FUNC_REL;
        optimizer_type = OptimizerType::LS;
        // set default parameters

        return ;
    }

    void OptimizerBase::update_params(Options *opt_in)
    {
        // update the parameters
    // max_iterations = opt_in -> max_iterations;
    // min_iterations = opt_in -> min_iterations;

    
        return ;
    }
}