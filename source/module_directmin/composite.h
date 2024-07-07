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

        Stiefel p1; // part 1
        Occupation p2; // part 2

        Composite& operator=(const Composite& comp);
        Composite& operator+=(const Composite& comp);

    };

} // end namespace ModuleDirectMin


#endif // COMPOSITE_H