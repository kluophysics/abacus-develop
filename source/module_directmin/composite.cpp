#include "composite.h"

namespace ModuleDirectMin
{

    Composite::Composite()
    {
        p1 = Stiefel(); // empty
        p2 = Occupation(); // empty
    }

    Composite::Composite(const Stiefel & s, const Occupation & occ)
    {
        p1 = s;
        p2 = occ;
        // return *this;
    }


    Composite & Composite::operator = (const Composite & comp)
    {
        p1 = comp.p1;
        p2 = comp.p2;
        return *this;

    }
    Composite & Composite::operator += (const Composite & comp)
    {
        p1 += comp.p1;
        p2 += comp.p2;
        return *this;

    }

    Composite Composite::t() 
    {
        p1 = this->p1.t();
        // p2 = this->p2;
        return *this;
    }

    // Composite Composite::operator*(const Composite & p) const
    // {
       
    // }

}