#include "optimizer_ls_rsd.h"


namespace Module_Optimizer
{
    
    RSD::RSD(Problem *prob_in, LSOptions *opt_in)
    {
        set_default_params();
        // OptimizerLSBase::initialize(prob_in);
        OptimizerLSBase::initialize();
        update_params(opt_in);
    }

    void RSD::set_default_params()
    {
        OptimizerLSBase::set_default_params();
        name.assign("RSD"); 

    }

    void RSD::update_params(LSOptions* opt_in)
    {
        OptimizerLSBase::update_params(opt_in);
        // set_cg_algo(opt_in->ls_cg_algo); // set the cg algorithm using input option
    }

    void RSD::get_search_direction()
    {
        d1 = -1.0 * gf1;
    }

    void RSD::update_data()
    {
        ;
    }


}