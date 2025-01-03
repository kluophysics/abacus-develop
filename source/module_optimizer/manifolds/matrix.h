#ifndef MATRIX_H
#define MATRIX_H


// define the product manifold of Matrix
#include <vector>
#include <armadillo>

namespace Module_Optimizer
{
    template<typename T>
    class Matrix;

    // template<typename T>
    // bool equal_dimension(const Matrix<T>& a, const Matrix<T>& b);

	//a Matrix class from armadillo
    template<typename T>
    class Matrix
    {
        public:
            Matrix();
            Matrix(const Matrix<T>&  var);
            Matrix(int nr, int nc);


            // arma::cx_mat operator[](const int k) const;
            // arma::cx_mat & operator[](const int k);


            // // uni operator overload
            // Matrix<T> operator-() const; // negative operator
            
            // // assignment operator overload
            // Matrix<T>& operator=(const Matrix<T>& p) ;
            // Matrix<T>& operator-=(const Matrix<T>& p);
            // Matrix<T>& operator+=(const Matrix<T>& p);
            // Matrix<T>& operator*=(const Matrix<T>& p);




            // Matrix<T> operator+( double s) const; // add p with s, p+s;
            // Matrix<T> operator-( double s) const; // p-s;
            // Matrix<T> operator*( double s) const; //  p*s;

            // Matrix<T> operator+( const Matrix<T> & p) const; // add p with s, p+s;
            // Matrix<T> operator-( const Matrix<T> & p) const; // add p with s, p+s;
            // Matrix<T> operator*( const Matrix<T> & p) const; // add p with s, p+s;

            // friend Matrix<T> operator+(double s, const Matrix<T> & p); // s+p
            // friend Matrix<T> operator*(double s, const Matrix<T> & p); // s*p
            // friend Matrix<T> operator-(double s, const Matrix<T> & p); // s-p

            
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

    // template<typename T>
    // bool equal_dimension(const Matrix<T> & a, const Matrix<T>  & b)
    // {
    //     return (a.nr == b.nr) && (a.nc == b.nc);
    // }
    // Friend function definition

}



#endif // MATRIX_H