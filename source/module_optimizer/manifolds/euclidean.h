#ifndef EUCLIDEAN_MANIFOLD_H
#define EUCLIDEAN_MANIFOLD_H

#include "manifold.h"
#include <armadillo>

namespace Module_Optimizer
{
    template<typename T>
    class Euclidean : public Manifold<T>
    {
    public:
        using ManifoldPoint = typename Manifold<T>::ManifoldPoint;
        using ManifoldVector = typename Manifold<T>::ManifoldVector;

        Euclidean(int nr, int nc, int num_manifolds) : 
            p(nr), n(nc), k(num_manifolds) {};
       

    private:
        int p; // Number of rows
        int n; // Number of columns
        int k; // Number of Stiefel manifolds in the product
    };
}

#endif // EUCLIDEAN_MANIFOLD_H