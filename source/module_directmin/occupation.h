#ifndef OCCUPATON_H
#define OCCUPATON_H

#include <vector>


namespace ModuleDirectMin
{
    class Occupation 
    {
        public:

        Occupation();
        Occupation(double * v, int l); // initialize from a double array

        ~Occupation();
    

        int size; // the size of the vector
        std::vector<double> occ_vector; // the vector of occ

        bool is_valid(); // check if the vector is valid

        
    
    };

} // namespace ModuleDirectMin


#endif // OCCUPATON_H