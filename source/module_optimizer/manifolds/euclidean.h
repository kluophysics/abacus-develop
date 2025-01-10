#ifndef EUCLIDEAN_MANIFOLD_H
#define EUCLIDEAN_MANIFOLD_H

#include "manifold.h"
#include <armadillo>

namespace Module_Optimizer
{
    class Euclidean : public Manifold
    {
    public:
        using ManifoldPoint = typename Manifold::ManifoldPoint;
        using ManifoldVector = typename Manifold::ManifoldVector;

        Euclidean(int nr, int nc, int num_manifolds) : 
            p(nr), n(nc), k(num_manifolds) {};
       

    private:
        int p; // Number of rows
        int n; // Number of columns
        int k; // Number of Stiefel manifolds in the product
    };
}

#endif // EUCLIDEAN_MANIFOLD_H