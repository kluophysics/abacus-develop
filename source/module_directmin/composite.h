#ifndef COMPOSITE_H
#define COMPOSITE_H

#include "stiefel.h"
#include "occupation.h"

// This header is to define the composite variable that 
// includes both the Stiefel variable and the occupation vector 
// By properly switching a flag, the occupation vector can be 
// unchanged, or update with the Stiefel variable (which is the case 
// in the finite temperature smearning )
// The design is not perfect but I hope this could serve as a 
// good starting point for future.

namespace ModuleDirectMin
{

    class Composite
    {
    public:

        Composite();
        // ~Composite();
        Composite(const Stiefel & s, const Occupation & occ);

        Stiefel p1; // part 1, Stiefel variable
        Occupation p2; // part 2, Occupation variable

        Composite& operator=(const Composite& comp);
        Composite& operator+=(const Composite& comp);
        Composite & operator*(const double s); // multiplicative by a factor s

        Composite t(); // complex transpose

        
        // Composite operator *(const Composite &p) const; // multiplicative operator
        
        // Occupation vector can change explicitly if true
        // fixed or only implicitly varied with the Stiefel variable if false
        bool explicit_occupation_flag;
        
        // apply a scaling factor beta to the direction vector for 
        // occupation variable, which might be helpful in convergence.
        Composite apply_async_scaling(double beta);
    };

} // end namespace ModuleDirectMin


#endif // COMPOSITE_H