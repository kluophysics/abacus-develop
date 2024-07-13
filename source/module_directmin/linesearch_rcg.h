#ifndef LINESEARCH_RCG_H
#define LINESEARCH_RCG_H

#include "linesearch_base.h"

namespace ModuleDirectMin
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


    class LineSearchRCG : public LineSearchBase
    {
    public:
        std::string *cg_algo_names;
        double sigma;

        CGAlgorithm cg_algo;
        ;

    protected:

        // update_data

        virtual void update_data();

        // set default params for line search
        virtual void set_default_parameters();

        // update params for line search using LineSearchOptions *ls_opt_in
        virtual void update_parameters(LineSearchOptions * ls_opt_in);
        
        //Compute the search direction based on the RCG forumla.
		//Reset the search direction to be negative gradient if the search direction is not sufficiently descent or
        //| g( grad f(x_{k+1}), grad f(x_{k}) ) | / g( grad f(x_{k+1}), grad f(x_{k+1}) ) > 0.1 (see [(5.52), NW06]]).
        virtual void get_search_direction() = 0;

        virtual void set_cg_algo( std::string& name);
        virtual void compute_sigma(); // compute sigma according to cg_algo;

    };
}

#endif /* LINESEARCH_RCG_H */