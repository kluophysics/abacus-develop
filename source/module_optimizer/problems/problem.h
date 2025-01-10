#ifndef PROBLEM_H
#define PROBLEM_H

// #include <armadillo>
#include <vector>
// #include "module_base/complexmatrix.h"
#include "manifolds/manifold.h"
namespace Module_Optimizer
{

    // class Manifold;

    // class Problem; when defining the problem, one has to give initial X0 and direction d0
    // template <typename T>
    class Problem
    {
    public:
        

        using ManifoldPoint = typename Manifold::ManifoldPoint;
        using ManifoldVector = typename Manifold::ManifoldVector;

        // Problem();
        virtual ~Problem() =0;

        // function value evaluated on x, which has to be defined for each problem
        // virtual double objective(const ManifoldPoint & x) = 0;
        virtual double obj(const ManifoldPoint & x)  const=0;

        // Euclidean gradient of f, which has to be defined for each problem
        // virtual ManifoldPoint grad(const ManifoldPoint& x) =0;
        virtual ManifoldVector & grad(const ManifoldPoint & x) const=0;

        // evaluate objective function and gradient at the same time.
        virtual void evaluate_obj_and_grad(const ManifoldPoint & x) const=0 ;

        // Riemannian Gradient defined on the manifold 
        virtual ManifoldVector & rie_grad(const ManifoldPoint & x ) const=0 ;

        // set the manifold of the objective function
        virtual void set_manifold(Manifold * mani_in) { mani = mani_in; }

        // get the manifold of the objective function
        inline Manifold * get_manifold() const { return mani; }  


        
        // variable
        Manifold * mani; // pointer to hold the manifold of the objective function 

    };

};

#endif