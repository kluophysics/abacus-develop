#ifndef SMARTSPACE_H
#define SMARTSPACE_H

#include <vector> // this is the core data structure use for Optim 

namespace Module_Optim
{
    // a vector of matrices denoted as SmartSpace class
    class SmartSpace
    {
        public:
            SmartSpace();
            SmartSpace(const SmartSpace& );
            /* SmartSpace supports up to three dimensions
            Default: it is a one-dimensional vector with n2=1, n3 = 1
            Else if: only 2 numbers are provided, then n3=1;
            Else: 3 numbers, it is a product of manifolds
            */
            SmartSpace(int n1, int n2=1, int n3=1);

        
            arma::cx_mat operator[](const int n3) const;
            arma::cx_mat & operator[](const int n3);


            // uni operator overload
            SmartSpace operator-() const; // negative operator
            
            // assignment operator overload
            SmartSpace& operator=(const SmartSpace& p) ;
            SmartSpace& operator-=(const SmartSpace& p);
            SmartSpace& operator+=(const SmartSpace& p);
            SmartSpace& operator*=(const SmartSpace& p);




            SmartSpace operator+( double s) const; // add p with s, p+s;
            SmartSpace operator-( double s) const; // p-s;
            SmartSpace operator*( double s) const; //  p*s;

            SmartSpace operator+( const SmartSpace & p) const; // add p with s, p+s;
            SmartSpace operator-( const SmartSpace & p) const; // add p with s, p+s;
            SmartSpace operator*( const SmartSpace & p) const; // add p with s, p+s;

            friend SmartSpace operator+(double s, const SmartSpace & p); // s+p
            friend SmartSpace operator*(double s, const SmartSpace & p); // s*p
            friend SmartSpace operator-(double s, const SmartSpace & p); // s-p

            
            void resize(int n3, int r, int c);

            double norm();
            SmartSpace t() const; // conjugate transpose.
            SmartSpace inv() const; // inverse.
            void print(const char* s) const;
            void brief_print(const char* s) const;

            void reset();
            void clear();

            inline int get_nn3() const { return nn3;}; 
            inline int get_nr() const { return nr;}; 
            inline int get_nc() const { return nc;}; 
            inline int get_size() const { return size;}; 

            // checn3 if two matrices a and b are of the same dimensions
            friend bool equal_dimension(const SmartSpace & a, const SmartSpace & b); // 


        private:
            int nn3; // number of manifolds, each manifold in this case is of same size
            int nr; // number of rows for each manifold
            int nc; // number of columns for each manifold
            int size; // nn3*nr*nc
            std::vector<arma::cx_mat> psm; // product of stiefel manifold
    };

}
#endif // SMARTSPACE_H
