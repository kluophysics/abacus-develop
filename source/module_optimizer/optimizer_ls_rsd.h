#ifndef OPTIMIZER_LS_RSD_H
#define OPTIMIZER_LS_RSD_H

#include "optimizer_ls_base.h"

namespace Module_Optimizer
{

    class RSD : public OptimizerLSBase
    {

    public:
        
        RSD(Problem *prob_in, LSOptions *opt_in);

        virtual void set_default_params();
        virtual void update_params(LSOptions* opt_in);

        virtual void get_search_direction();
        virtual void update_data();

    };

}

#endif // OPTIMIZER_LS_RSD_H