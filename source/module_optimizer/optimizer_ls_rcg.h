#ifndef  OPTIMIZER_LS_RCG_H
#define OPTIMIZER_LS_RCG_H

#include "optimizer_ls_base.h"

namespace Module_Optimizer
{

    enum CGAlgorithm {
        FLETCHER_REEVES,
        DAI_YUAN,
        POLAK_RIBIERES,
        POLAK_RIBIERES_MOD,
        FR_PR,
        HAGER_ZHANG,
        HESTENES_STIEFEL,
        RCG_ALGO_LENGTH
    };

    class RCG : public OptimizerLSBase
    {
    public:
        // RCG(); // default constructor, no specified anything
        // RCG(Problem *prob);

        // constructor
        RCG(Problem *prob_in, LSOptions * opt_in);
        // ~RCG();

        CGAlgorithm cg_algo;
        LSOptions *opt; // options used

    protected:
    
        // void optimize();
        virtual void update_data();
        virtual void set_default_params();
        virtual void update_params(LSOptions* opt_in);

        //Compute the search direction based on the RCG forumla.
		//Reset the search direction to be negative gradient if the search direction is not sufficiently descent or
        //| g( grad f(x_{k+1}), grad f(x_{k}) ) | / g( grad f(x_{k+1}), grad f(x_{k+1}) ) > 0.1 (see [(5.52), NW06]]).
        virtual void get_search_direction();

        virtual void set_cg_algo( std::string& name);
        virtual void compute_sigma(); // compute sigma according to cg_algo;

        // virtual void print_info(); // print info related RCG

        std::string *cg_algo_names;

        double sigma;
        /*sigma is the coefficient in - \grad f(x_{k+1}) + \sigma \eta_k*/
    };


} // end namespace Module_Optimizer


#endif // OPTIMIZER_LS_RCG_H