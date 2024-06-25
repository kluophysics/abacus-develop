#ifndef STIEFEL_H
#define STIEFEL_H

// #include "module_base/complexmatrix.h"
// define the product manifold of Stiefel
#include <vector>
#include <armadillo>

namespace ModuleDirectMin
{

    class Stiefel
    {
        public:
            Stiefel();
            Stiefel(const Stiefel& );
            Stiefel(int k, int r, int c);

            // Stiefel operator+=(const Stiefel&p);

            // Stiefel operator*( double s);

            // assignment operator overload
            Stiefel& operator=(const Stiefel& p);
            Stiefel& operator-=(const Stiefel& p);
            Stiefel& operator+=(const Stiefel& p);
            Stiefel& operator*=(const Stiefel& p);

            // ModuleBase::ComplexMatrix &operator()(int k);
            arma::cx_mat operator[](const int k) const;
            arma::cx_mat & operator[](const int k);

            void resize(int k, int r, int c);

            double norm();
            Stiefel t() const; // conjugate transpose.
            Stiefel inv() const; // inverse.
            void print(const char* s) const;
            void brief_print(const char* s) const;

            void reset();
            void clear();

        // private:
            int nk; // number of manifolds, each manifold in this case is of same size

            int nr; // number of columns for each manifold
            int nc; // number of columns for each manifold
            
            int size; // nk*nr*nc
            // std::vector<ModuleBase::ComplexMatrix> psm;
            std::vector<arma::cx_mat> psm;

            
    };

    // uni operator overload
    const Stiefel operator-(const Stiefel&p); // negative operator
    // Stiefel  operator+();
    
    // binary operator overload
    const Stiefel operator+(const Stiefel&p,const   Stiefel& q); // addition operator
    const Stiefel operator-(const Stiefel&p,const   Stiefel& q); // addition operator
    Stiefel operator*(const Stiefel&p,const   Stiefel& q); // multiplicative operator
    Stiefel operator%(const Stiefel&p,const   Stiefel& q); // elementwise
    Stiefel operator/(const Stiefel&p,const   Stiefel& q); // elementwise division

    Stiefel operator+( double s, const Stiefel&p); // add p with s
    Stiefel operator+( const Stiefel&p,  double s); // add p with s
    Stiefel operator-( double s, const Stiefel&p); //  p with s
    Stiefel operator-( const Stiefel&p,  double s); // add p with s

    Stiefel operator*( double s, const Stiefel&p); // scale p with s
    Stiefel operator*(const Stiefel&p,  double s); // same as above but order is reversed.

    // Projection onto the tangent space to the Stiefel manifold.
    Stiefel proj(const Stiefel& X, const Stiefel& Z) ;

    // vector transport in Stiefel manifold
    Stiefel vectran(const Stiefel& X, const Stiefel& Z);

    // retraction function in Stiefel manifold
    Stiefel retraction(const Stiefel& X, const Stiefel& Z);

    // differentiated retraction function in Stiefel manifold
    Stiefel diff_retraction(const Stiefel& X, const Stiefel& Z);

    // canonical metric for stiefel manifold
    double metric(const Stiefel& X, const Stiefel& A, const Stiefel& B);

    bool is_orthogonal(const Stiefel& X);
}


#endif