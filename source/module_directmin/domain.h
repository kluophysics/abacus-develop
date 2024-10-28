#ifndef DOMAIN_H
#define DOMAIN_H

#define DomainVector DomainPoint

#include "manifold/stiefel.h"
#include "occupation.h"

namespace ModuleDirectMin
{
    class DomainPoint;
    // class DomainVector;

    class Domain : public Stiefel
    {
    public:
        // Occupation vector can change explicitly if true
        // fixed or only implicitly varied with the Stiefel variable if false
        bool explicit_occupation_flag;

        // apply a scaling factor beta to the direction vector for 
        // occupation variable, which might be helpful in convergence.
        void apply_async_scaling(DomainVector * dv, double beta);


        // virtual DomainPoint retraction(const DomainPoint& x, const DomainVector& d);
        // virtual DomainVector diff_retraction(const DomainPoint &x, const DomainVector d1,
        //                     const DomainPoint & x2, const DomainVector & d2);
        virtual double norm(const DomainPoint &x, 
                 DomainVector& d1, const DomainVector& d2);
    };

    class DomainPoint : public StiefelPoint
    {
    public:
        // DomainPoint(const DomainPoint & dp);
        DomainPoint(const DomainPoint & dp);
        DomainPoint(StiefelPoint & x_in, const Occupation & occ_in);


        StiefelPoint x;
        Occupation occ; // Occupation variable
        
        friend bool equal_dimension(const DomainPoint & a, const DomainPoint & b); // 
        // friend bool equal_dimension(const DomainPoint & a, const DomainVector & b); // 
        // friend bool equal_dimension(const DomainVector & a, const DomainVector & b); // 
        // friend bool equal_dimension(const DomainVector & a, const DomainVector & b); // 


    };

    // class DomainVector:
}



#endif // DOMAIN_H