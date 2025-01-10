#ifndef OPTIMIZER_LS_BASE_H
#define OPTIMIZER_LS_BASE_H

#include "manifold.h"

namespace Module_Optimizer {

    class Optimizer_LS_Base {
    public:
        void Armijo();
        void Wolfe();
        void StrongWolfe();

    private:
        ManifoldPoint current_point;
        ManifoldVector current_gradient;
        ManifoldVector search_direction;
        double step_size;

        double objective_function(const ManifoldPoint &x);
        ManifoldVector gradient(const ManifoldPoint &x);
        Manifold manifold;
    };

}

#endif // OPTIMIZER_LS_BASE_H