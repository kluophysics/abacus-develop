#ifndef VECTOR_H
#define VECTOR_H

#include <vector> // this is the core data structure use for Optim 

namespace Module_Optim
{
    // a vector of matrices denoted as OptimVector class
    class OptimVector
    {
        public:
            OptimVector();
            OptimVector(const OptimVector& );
            /* OptimVector supports up to three dimensions
            Default: it is a one-dimensional vector with n2=1, n3 = 1
            Else if: only 2 numbers are provided, then n3=1;
            Else: 3 numbers, it is a product of manifolds
            */
            OptimVector(int n1, int n2=1, int n3=1);

        
            arma::cx_mat operator[](const int n3) const;
            arma::cx_mat & operator[](const int n3);


            // uni operator overload
            OptimVector operator-() const; // negative operator
            
            // assignment operator overload
            OptimVector& operator=(const OptimVector& p) ;
            OptimVector& operator-=(const OptimVector& p);
            OptimVector& operator+=(const OptimVector& p);
            OptimVector& operator*=(const OptimVector& p);




            OptimVector operator+( double s) const; // add p with s, p+s;
            OptimVector operator-( double s) const; // p-s;
            OptimVector operator*( double s) const; //  p*s;

            OptimVector operator+( const OptimVector & p) const; // add p with s, p+s;
            OptimVector operator-( const OptimVector & p) const; // add p with s, p+s;
            OptimVector operator*( const OptimVector & p) const; // add p with s, p+s;

            friend OptimVector operator+(double s, const OptimVector & p); // s+p
            friend OptimVector operator*(double s, const OptimVector & p); // s*p
            friend OptimVector operator-(double s, const OptimVector & p); // s-p

            
            void resize(int n3, int r, int c);

            double norm();
            OptimVector t() const; // conjugate transpose.
            OptimVector inv() const; // inverse.
            void print(const char* s) const;
            void brief_print(const char* s) const;

            void reset();
            void clear();

            inline int get_nn3() const { return nn3;}; 
            inline int get_nr() const { return nr;}; 
            inline int get_nc() const { return nc;}; 
            inline int get_size() const { return size;}; 

            // checn3 if two matrices a and b are of the same dimensions
            friend bool equal_dimension(const OptimVector & a, const OptimVector & b); // 


        private:
            int nn3; // number of manifolds, each manifold in this case is of same size
            int nr; // number of rows for each manifold
            int nc; // number of columns for each manifold
            int size; // nn3*nr*nc
            std::vector<arma::cx_mat> psm; // product of stiefel manifold
    };

}
#endif // VECTOR_H
