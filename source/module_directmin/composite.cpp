#include "composite.h"
#include <cassert>

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
    Composite Composite::operator - () const
    {
        // p1 = - p1;
        // p2.occ_vector = - p2.occ_vector;
        return Composite(-p1, -p2);
    }
    Composite & Composite::operator*(double s)
    {
        p1 = this->p1 * s;
        p2 = this->p2 * s;
        return *this;
    }

    Composite & operator*(double s, const Composite & comp)
    {
        return s*comp;
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

    double Composite::inner_product(const Composite &c1, const Composite &c2)
    {
        double total_metric = 0.0;
        // check size consistency
        assert(this->p1.size == c1.p1.size);
        assert(this->p2.size == c1.p2.size);
        assert(c1.p1.size == c2.p1.size);
        assert(c1.p2.size == c2.p2.size);

        // metric contribution from the occupation vector
        if(explicit_occupation_flag)
        {
            for(int i=0; i < c1.p2.size; ++i)
            {
                total_metric += c1.p2.occ_vector[i] * c2.p2.occ_vector[i];
            }
        }

        // metric contribution from the Stiefel variable
        total_metric += (this->p1).metric(c1.p1, c2.p1);

        return total_metric;
    }

    // Composite Composite::operator*(const Composite & p) const
    // {
       
    // }
    Composite Composite::apply_async_scaling(double beta=1.0)
    {
        p2 = beta * this->p2;
        return * this;
    }

}