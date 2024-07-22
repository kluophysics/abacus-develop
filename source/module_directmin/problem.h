#ifndef PROBLEM_H
#define PROBLEM_H

#include "domain.h"

namespace ModuleDirectMin
{
    class Problem
    {
    public:
        
        // Problem();
        virtual ~Problem();

        // function value evaluated on x, which has to be defined for each problem
        virtual double fval(const Domain& x) = 0;
        // function value evaluated on matrix x and vector v, made for occupation as vector v;
        // virtual double fval(const Domain& x,  arma::vec & v) =0;


        // Euclidean gradient of fval, which has to be defined for each problem
        virtual Domain grad(const Domain& x) =0;

        // Riemannian gradient defined on the Stiefel manifold 
        virtual Domain rgrad(const Domain& x);

        // preconditioner...
        virtual Domain  preconditioner(const Domain & C_in);


        virtual void set_domain(Domain * domain_in);
        
        /*Obtain the domain manifold of the cost function*/
        inline Domain * get_domain() const { return domain; };
        // print information
        // virtual void print();
    protected:
        Domain * domain;
    };

}

#endif /* PROBLEM_H */