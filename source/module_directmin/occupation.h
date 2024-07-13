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

        Occupation& operator=(const Occupation& occ);
        Occupation& operator+=(const Occupation& occ);
        Occupation& operator-=(const Occupation& occ);
        // Occupation& operator*=(const Occupation& occ);
        // Occupation& operator/=(const Occupation& occ);

        Occupation operator - () const; // negation operator
        Occupation operator*( double s) const; // add p with s
        friend Occupation operator*(double s, const Occupation& occ);

        
    };



} // namespace ModuleDirectMin


#endif // OCCUPATON_H