#ifndef PROBLEM_H
#define PROBLEM_H

// #include <armadillo>
#include <vector>
// #include "module_base/complexmatrix.h"
#include "../manifolds/stiefel.h"
namespace Module_Optimizer
{

    // class Problem; when defining the problem, one has to give initial X0 and direction d0
    template <typename T>
    class Problem
    {
    public:
        
        using typename Manifold<T>::ManifoldPoint;
        using typename Manifold<T>::ManifoldVector;

        // Problem();
        virtual ~Problem();

        // function value evaluated on x, which has to be defined for each problem
        // virtual double objective(const ManifoldPoint & x) = 0;
        virtual double obj() ;

        // Euclidean gradient of f, which has to be defined for each problem
        // virtual ManifoldPoint grad(const ManifoldPoint& x) =0;
        virtual ManifoldPoint grad();

        // evaluate objective function and gradient at the same time.
        virtual void evaluate_obj_and_grad();

        // Riemannian Gradient defined on the manifold 
        virtual ManifoldVector RieGrad();

        // initial value for X 
        ManifoldPoint X0;
        // initial value for d
        ManifoldVector d0;
        
        // variable
        // ManifoldPoint X; 
    };

};

#endif