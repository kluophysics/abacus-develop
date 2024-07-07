#include "occupation.h"

namespace ModuleDirectMin
{

    Occupation::Occupation()
    {
        size = 0;
        occ_vector.resize(size);
    }

    Occupation::Occupation(double *v, int l)
    {
        size = l;
        occ_vector.resize(size);

        for(int i =0; i< size; i++)
        {
            occ_vector[i] = v[i];
        }
    }

    Occupation::~Occupation()
    {
        size = 0;
        occ_vector.clear();
    }


    bool Occupation::is_valid() 
    {
        for(int i = 0; i < size; i++)
        {
            if(occ_vector[i] < 0.0) // occupation has to be greater than zero
            {
                return false;
            } 
        }
        return true;
    }
} // namespace ModuleDirectMin