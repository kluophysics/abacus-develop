#include "occupation.h"
#include <cassert>

namespace ModuleDirectMin
{

    Occupation::Occupation()
    {
        norb = 1;
        nk = 1;
        size = nk * norb;
        occ_vector.resize(nk);
    }
    Occupation::Occupation(int nk_in, int norb_in)
    {
        norb = norb_in;
        nk = nk_in;
        size = nk * norb;
        for(int ik=0; ik < nk; ik++)
        {
            occ_vector[ik].resize(norb);
        }
    }
    // Occupation::Occupation(double *v, int l)
    // {
    //     size = l;
    //     occ_vector.resize(size);

    //     for(int i =0; i< size; i++)
    //     {
    //         occ_vector[i] = v[i];
    //     }
    // }

    Occupation::~Occupation()
    {
        // size = 0;
        // nk = 0;
        for(int ik=0; ik < this->nk; ik++)
        {
            occ_vector[ik].clear();
        }
        
    }


    bool Occupation::is_valid() 
    {
        for(int ik = 0; ik < nk; ik ++)
        {
            for(int i = 0; i < size; i++)
            {
                if((occ_vector[ik]) (i) < 0.0) // occupation has to be greater than zero
                {
                    return false;
                } 
            }
        }

        return true;
    }

    Occupation& Occupation::operator=(const Occupation& occ)
    {
        size = occ.size;
        nk = occ.nk;
        norb = occ.norb;
        occ_vector = occ.occ_vector;
        
        return *this;


    }

    Occupation& Occupation::operator+=(const Occupation& occ)
    {
        assert(this->size == occ.size);
        assert(this->nk == occ.nk);
        assert(this->norb = occ.norb);

        for(int ik=0; ik < this->nk; ik++)
        {
            occ_vector[ik] += occ.occ_vector[ik];
        }
        return *this;
    }
    Occupation& Occupation::operator-=(const Occupation& occ)
    {
        assert(this->size == occ.size);
        assert(this->nk == occ.nk);
        assert(this->norb = occ.norb);

        for(int ik=0; ik < this->nk; ik++)
        {
            occ_vector[ik] -= occ.occ_vector[ik];
        }
        return *this;
    }

    Occupation Occupation::operator-() const
    {
        Occupation result = *this;

        for(int ik=0; ik < this->nk; ik++)
        {
            result.occ_vector[ik] = -this->occ_vector[ik];
        }
        return result;
    }
    
    Occupation Occupation::operator*(double s) const
    {
        Occupation result = *this;

        for(int ik=0; ik < this->nk; ik++)
        {
            result.occ_vector[ik] = s*this->occ_vector[ik];
        }
        return result;
    }
    Occupation operator*(double s, const Occupation & occ)
    {
        Occupation result(occ);
        for(int ik=0; ik < occ.nk; ik++)
        {
            result.occ_vector[ik] = s*occ.occ_vector[ik];
        }
        return result;
    }


} // namespace ModuleDirectMin