#ifndef OCCUPATON_H
#define OCCUPATON_H

#include <vector>
#include <armadillo>


namespace ModuleDirectMin
{
    class Occupation 
    {
        public:

        Occupation();
        Occupation(int nk_in, int norb_in);
        Occupation(double ** v, int nk_in, int norb_in); // initialize from a double array

        ~Occupation();
    

        int size; // the total  size
        int norb; // the number of orbitals
        int nk; // the number of k points
        std::vector<arma::vec> occ_vector; // the vector of occ

        bool is_valid(); // check if the vector is valid

        Occupation& operator=(const Occupation& occ);
        Occupation& operator+=(const Occupation& occ);
        Occupation& operator-=(const Occupation& occ);
        // Occupation& operator*=(const Occupation& occ);
        // Occupation& operator/=(const Occupation& occ);

        Occupation operator - () const; // negation operator
        Occupation operator*( double s) const; // add p with s
        friend Occupation operator*(double s, const Occupation& occ);

        // double operator*(const Occupation& occ) const;
        // friend Occupation operator*(const Occupation& a, const Occupation & b) const;

        
    };



} // namespace ModuleDirectMin


#endif // OCCUPATON_H