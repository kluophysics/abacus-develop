#ifndef STIEFEL_H
#define STIEFEL_H

#include "options.h"

// #include "module_base/complexmatrix.h"
// define the product manifold of Stiefel
#include <vector>
#include <armadillo>

namespace ModuleDirectMin
{
	/*Note that not all metrics, retractions and vector transports have been done.*/

	// /* Riemannian Metric for the Stiefel manifold:
	// Eucldean: g_x(etax, xix) = \trace(etax^T xix);
	// Canonical: g_x(etax, xix) = \trace(etax^T (I_n - x x^T / 2) xix); */
	// enum MetricType{ STIE_EUCLIDEAN, STIE_CANONICAL, STIEMETRICLENGTH };

	// /*Retraction for the Stiefel manifold
	// RT_QF: qf retraction defined in [AMS2008, (4.8)]
	// RT_POLAR: polar based retraction defined in [AMS2008, (4.7)]
	// RT_EXP: The exponential mapping
	// CAYLEYR: the Cayley transform in [Zhu2016]
	// [AMS2008]P.-A. Absil, R. Mahony, and R. Sepulchre. Optimization algorithms on matrix manifolds.
	// Princeton University Press, Princeton, NJ, 2008.
	// [HGA2015]:Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
	// SIAM Journal on Optimization, 25(3):1660?685,2015.
	// [Hua2013]:W. Huang. Optimization algorithms on Riemannian manifolds with applications.
	// PhD thesis, Florida State University, Department of Mathematics, 2013.
	// [Zhu2016]: Xiaojing Zhu, A Riemannian conjugate gradient method for optimization on the Stiefel Manifold */
	// enum RetractionType{ STIE_QF, STIE_POLAR, STIE_EXP, STIE_CAYLEYR, STIERETRACTIONLENGTH };

	// /*Vector transport for the Stiefel manifold
	// VT_PARALLELIZATION: Vector transport by parallelization, See [HAG2015, Section 2.3.1]
	// VT_RIGGING: Vector transport by rigging, See [HAG2015, Section 2.3.2]
	// VT_PARALLELTRANSLATION: parallel translation
	// VT_CAYLEY: the vector transport based on Cayley transform. [Zhu2016]
	// [HAG2015]:W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trust-region method.
	// Mathematical Programming, 150(2):179?16, February 2015
	// [Zhu2016]: Xiaojing Zhu, A Riemannian conjugate gradient method for optimization on the Stiefel Manifold */
	// enum VectorTranportType{ STIE_PARALLELIZATION, STIE_RIGGING, STIE_PARALLELTRANSLATION, STIE_PROJECTION, STIE_CAYLEYVT, STIEVECTORTRANSPORTLENGTH };

    class Stiefel
    {
        public:
            Stiefel();
            Stiefel(const Stiefel& );
            Stiefel(int k, int r, int c);

            // Stiefel operator+=(const Stiefel&p);

            // Stiefel operator*( double s);
            // uni operator overload
            Stiefel operator-() const; // negative operator
            
            // assignment operator overload
            Stiefel& operator=(const Stiefel& p);
            Stiefel& operator-=(const Stiefel& p);
            Stiefel& operator+=(const Stiefel& p);
            Stiefel& operator*=(const Stiefel& p);

            // ModuleBase::ComplexMatrix &operator()(int k);
            arma::cx_mat operator[](const int k) const;
            arma::cx_mat & operator[](const int k);


            Stiefel operator+( double s) const; // add p with s, p+s;
            Stiefel operator-( double s) const; // p-s;
            Stiefel operator*( double s) const; //  p*s;

            Stiefel operator+( const Stiefel & p) const; // add p with s, p+s;
            Stiefel operator-( const Stiefel & p) const; // add p with s, p+s;
            Stiefel operator*( const Stiefel & p) const; // add p with s, p+s;

            friend Stiefel operator+(double s, const Stiefel & p); // s+p
            friend Stiefel operator*(double s, const Stiefel & p); // s*p
            friend Stiefel operator-(double s, const Stiefel & p); // s-p

            
            void resize(int k, int r, int c);

            double norm();
            Stiefel t() const; // conjugate transpose.
            Stiefel inv() const; // inverse.
            void print(const char* s) const;
            void brief_print(const char* s) const;

            void reset();
            void clear();

            int nk; // number of manifolds, each manifold in this case is of same size

            int nr; // number of rows for each manifold
            int nc; // number of columns for each manifold
            
            int size; // nk*nr*nc
            // std::vector<ModuleBase::ComplexMatrix> psm;
            std::vector<arma::cx_mat> psm; // product of stiefel manifold


        // protected:
            MetricType metric_type; // Riemannian metric
            RetractionType retraction_type; // retraction method
            VectorTranportType vector_transport_type; // vector transport method

            void switch_metric(MetricType metric);
            void switch_retraction(RetractionType retraction_type);
            void switch_vector_transport(VectorTranportType vector_transport_type);


            // Stiefel projection(const Stiefel& Z);     // Projection onto the tangent space to the Stiefel manifold.
            // Stiefel vector_transport(const Stiefel& Z);     // vector transport in Stiefel manifold
            // Stiefel retraction(const Stiefel& Z);             // retraction function in Stiefel manifold
            
            // Stiefel diff_retraction(const Stiefel& Z);     // differentiated retraction function in Stiefel manifold
            
            double metric(const Stiefel& A, const Stiefel& B); // metric for stiefel manifold
            double canonical_metric(const Stiefel& A, const Stiefel& B); // canonical metric for stiefel manifold
            double euclidean_metric(const Stiefel& A, const Stiefel& B); // euclidean metric
            bool is_orthogonal(); // check for orthogonal

            // friend class Options; // methods allowed to be changed from Options

    };


    // Stiefel  operator+();
    
    // binary operator overload
    // const Stiefel operator+(const Stiefel&p,const   Stiefel& q); // addition operator
    // const Stiefel operator-(const Stiefel&p,const   Stiefel& q); // addition operator
    // Stiefel operator*(const Stiefel&p,const   Stiefel& q); // multiplicative operator
    // Stiefel operator%(const Stiefel&p,const   Stiefel& q); // elementwise
    // Stiefel operator/(const Stiefel&p,const   Stiefel& q); // elementwise division

    // Stiefel operator+( double s, const Stiefel&p); // add p with s
    // Stiefel operator+( const Stiefel&p,  double s); // add p with s
    // Stiefel operator-( double s, const Stiefel&p); //  p with s
    // Stiefel operator-( const Stiefel&p,  double s); // add p with s

    // Stiefel operator*( double s, const Stiefel&p); // scale p with s
    // Stiefel operator*(const Stiefel&p,  double s); // same as above but order is reversed.

    // Projection onto the tangent space to the Stiefel manifold.
    // Stiefel proj(const Stiefel& X, const Stiefel& Z) ;

    // vector transport in Stiefel manifold
    Stiefel vector_transport(const Stiefel& x, const Stiefel& etax,
    const Stiefel& y, const Stiefel& xix);

    // // retraction function in Stiefel manifold
    Stiefel retraction(const Stiefel& x, const Stiefel& etax);

    // // differentiated retraction function in Stiefel manifold
    // Stiefel diff_retraction(const Stiefel& X, const Stiefel& Z);

    // canonical metric for stiefel manifold
    double metric(const Stiefel& X, const Stiefel& A, const Stiefel& B);

    // bool is_orthogonal(const Stiefel& X);
}


#endif