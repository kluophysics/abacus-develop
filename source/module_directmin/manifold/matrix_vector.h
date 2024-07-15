#ifndef MATRIX_VECTOR_H
#define MATRIX_VECTOR_H


// define the product manifold of MatrixVector
#include <vector>
#include <armadillo>

namespace ModuleDirectMin
{
	
    class MatrixVector
    {
        public:
            MatrixVector();
            MatrixVector(const MatrixVector& );
            MatrixVector(int k, int r, int c);


            arma::cx_mat operator[](const int k) const;
            arma::cx_mat & operator[](const int k);


            // uni operator overload
            MatrixVector operator-() const; // negative operator
            
            // assignment operator overload
            MatrixVector& operator=(const MatrixVector& p);
            MatrixVector& operator-=(const MatrixVector& p);
            MatrixVector& operator+=(const MatrixVector& p);
            MatrixVector& operator*=(const MatrixVector& p);




            MatrixVector operator+( double s) const; // add p with s, p+s;
            MatrixVector operator-( double s) const; // p-s;
            MatrixVector operator*( double s) const; //  p*s;

            MatrixVector operator+( const MatrixVector & p) const; // add p with s, p+s;
            MatrixVector operator-( const MatrixVector & p) const; // add p with s, p+s;
            MatrixVector operator*( const MatrixVector & p) const; // add p with s, p+s;

            friend MatrixVector operator+(double s, const MatrixVector & p); // s+p
            friend MatrixVector operator*(double s, const MatrixVector & p); // s*p
            friend MatrixVector operator-(double s, const MatrixVector & p); // s-p

            
            void resize(int k, int r, int c);

            double norm();
            MatrixVector t() const; // conjugate transpose.
            MatrixVector inv() const; // inverse.
            void print(const char* s) const;
            void brief_print(const char* s) const;

            void reset();
            void clear();

            inline int get_nk() { return nk;}; 
            inline int get_nr() { return nr;}; 
            inline int get_nc() { return nc;}; 
            inline int get_size() { return size;}; 



        private:
            int nk; // number of manifolds, each manifold in this case is of same size
            int nr; // number of rows for each manifold
            int nc; // number of columns for each manifold
            int size; // nk*nr*nc
            std::vector<arma::cx_mat> psm; // product of stiefel manifold
    };
}



#endif // MATRIX_VECTOR_H