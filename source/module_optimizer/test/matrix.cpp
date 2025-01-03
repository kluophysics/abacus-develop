#include "matrix.h"
#include <string>
#include <cassert>
#include <complex>

namespace Module_Optimizer
{

    template<typename T>
	Matrix<T>::Matrix()
    {
        nr = 0;
        nc = 0;
        size = 0;
        arma::Mat<T> data(nr, nc);
    }

    template<typename T>
	Matrix<T>::Matrix(const Matrix& var)
    {
        nr = var.nr;
        nc = var.nc;
        size = var.size;
        arma::Mat<T> data(nr, nc);

    }

    template<typename T>
	Matrix<T>::Matrix( int nr_, int nc_)
    {
        nr = nr_;
        nc = nc_;
        size = nr*nc;
        arma::Mat<T> data(nr, nc);
    }

    // template class Matrix<float>;
    template class Matrix<double>;
    template class Matrix<std::complex<double>>;

};
