#ifndef PROBLEM_H
#define PROBLEM_H

#include "composite.h"

namespace ModuleDirectMin
{
    class Problem
    {
    public:
        
        // Problem();
        virtual ~Problem();

        // function value evaluated on x, which has to be defined for each problem
        virtual double fval(const Composite& x) = 0;
        // function value evaluated on matrix x and vector v, made for occupation as vector v;
        // virtual double fval(const Composite& x,  arma::vec & v) =0;


        // Euclidean gradient of fval, which has to be defined for each problem
        virtual Composite grad(const Composite& x) =0;

        // Riemannian gradient defined on the Stiefel manifold 
        virtual Composite rgrad(const Composite& x);

        // preconditioner...
        virtual Composite  preconditioner(const Composite & C_in);


        virtual void set_domain(Composite * domain_in);
        
        // print information
        // virtual void print();
    protected:
        Composite * domain;
    };

}

#endif /* PROBLEM_H */