#ifndef MATRIX_H
#define MATRIX_H


// define the product manifold of Matrix
#include <vector>
#include <armadillo>

namespace Module_Optimizer
{
	//a Matrix class from armadillo
    template<typename T>
    class Matrix
    {
        public:
            Matrix();
            Matrix(const Matrix<T>&  var);
            Matrix(int nr, int nc);

            // void resize(int nr, int nc);

            // double norm();
            // Matrix<T> t() const; // conjugate transpose.
            // Matrix<T> inv() const; // inverse.
            // void print(const char* s) const;
            // void brief_print(const char* s) const;

            // void reset();
            // void clear();

            inline int get_nr() const { return nr;}; 
            inline int get_nc() const { return nc;}; 
            inline int get_size() const { return size;}; 
            inline  double norm() const { return arma::norm(data);}; // norm.
            inline Matrix<T> inverse() const {return arma::inv(data); } ;// inverse.
            inline Matrix<T> conjugate() const {return arma::conj(data); } ;//  conjugate transpose.
            inline Matrix<T> transpose() const {return arma::trans(data); } ;//  complex conjugate and  transpose.

            //  friend bool equal_dimension<>(const Matrix<T>& a, const Matrix<T>& b);


        private:
            int nr; // number of rows for each manifold
            int nc; // number of columns for each manifold
            int size; // nr*nc

            bool is_complex; // check if the matrix is complex

            arma::Mat<T> data; // data of the matrix
    };

}


#endif // MATRIX_H