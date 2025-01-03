#include "matrix.h"
#include <string>
#include <cassert>

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


    // template<typename T>
    // bool equal_dimension(const Matrix<T>& a, const Matrix<T>& b) {
    //     return (a.nr == b.nr) && (a.nc == b.nc);
    // }

    //  template<typename T>
	// Matrix<T>::resize(int k, int r, int c)
    // {
    //     psm.resize(k);

    //     size = k*r*c;
    //     nk = k;
    //     nr = r;
    //     nc = c;
    //     for(int ik = 0; ik < k; ik++)
    //     {
    //         psm[ik] = arma::cx_mat(nr, nc, arma::fill::zeros);
    //     }
    // }

    // void template<typename T>
	// Matrix<T>::reset()
    // {
    //     for( int ik = 0; ik < nk; ik++)
    //     {
    //         psm[ik].zeros();
    //     }
    // }

    // void template<typename T>
	// Matrix<T>::clear()
    // {
    //     nk = 0;
    //     nr = 0;
    //     nc = 0;
    //     size = 0;
    //     for( int ik = 0; ik < nk; ik++)
    //     {
    //         psm[ik].clear();
    //     }
    // }


    // arma::cx_mat template<typename T>
	// Matrix<T>::operator[](const int k) const
    // {
    //     assert( k >= 0 && k < nk);
    //     return psm[k];
    // }

    // arma::cx_mat & template<typename T>
	// Matrix<T>::operator[](const int k) 
    // {
    //     assert( k >= 0 && k < nk);
    //     return psm[k];
    // }


    // template<typename T>
	//  Matrix<T> & Matrix<T>::operator-() const
    // {
    //     // Matrix result(p);

    //     for(int ik = 0; ik < this->nk; ik++)
    //     {
    //         (*this)[ik] = -(*this)[ik] ;
    //     }
    //     return *this;

    // }


    // template<typename T>
	// Matrix<T>& Matrix<T>::operator=(const Matrix&p)
    // {
    //     this->resize(p.nk, p.nr, p.nc);

    //     for(int ik = 0; ik < p.nk; ik++)
    //     {
    //         (*this)[ik] = p[ik];
    //     }
    //     return *this;
    // }

    // template<typename T>
	//  Matrix<T>& Matrix<T>::operator-=(const Matrix&p)
    // {
    //     this->resize(p.nk, p.nr, p.nc);

    //     for(int ik = 0; ik < p.nk; ik++)
    //     {
    //         (*this)[ik] -= p[ik];
    //     }
    //     return *this;
    // }

    // template<typename T>
	//  Matrix<T>& Matrix<T>::operator+=(const Matrix&p)
    // {
    //     this->resize(p.nk, p.nr, p.nc);

    //     for(int ik = 0; ik < p.nk; ik++)
    //     {
    //         (*this)[ik] += p[ik];
    //     }
    //     return *this;
    // }

//     template<typename T> 
//     Matrix<T>& Matrix<T>::operator*=(const Matrix&p)
//     {
//         assert(this->nk == p.nk);
//         assert(this->nc == p.nr);

//         this->resize(p.nk, this->nr, p.nc);

//         for(int ik = 0; ik < p.nk; ik++)
//         {
//             (*this)[ik] *= p[ik];
//         }
//         return *this;
//     }

//    template<typename T>
// 	 Matrix & Matrix<T>::operator +(double s) const
//     {
//         for(int ik = 0; ik < this->nk; ik++)
//         {
//             (*this)[ik] += s;
//         }
//         return *this;
//     }

//     template<typename T>
// 	Matrix & Matrix<T>::operator *(double s) const
//     {
//         for(int ik = 0; ik < this->nk; ik++)
//         {
//             (*this)[ik] *= s;
//         }
//         return *this;
//     }

//    template<typename T>
// 	 Matrix<T> Matrix<T>::operator -(double s) const
//     {
//         for(int ik = 0; ik < this->nk; ik++)
//         {
//             (*this)[ik] -= s;
//         }
//         return *this;
//     }

//     template<typename T>
//     Matrix  <T> Matrix<T>::operator+( const Matrix & p) const
//     {
//         assert(this->nk == p.nk);
//         assert(this->nc == p.nc);
//         assert(this->nr == p.nr);

//         for(int ik = 0; ik < this->nk; ik++)
//         {
//             (*this)[ik] += p[ik];
//         }
//         return *this;
//     }
//     template<typename T>
// 	Matrix<T> Matrix<T>::operator-( const Matrix & p) const
//     {
//         assert(this->nk == p.nk);
//         assert(this->nc == p.nc);
//         assert(this->nr == p.nr);

//         for(int ik = 0; ik < this->nk; ik++)
//         {
//             (*this)[ik] -= p[ik];
//         }
//         return *this;
//     }

//     template<typename T>
// 	Matrix<T> Matrix<T>::operator*( const Matrix & p) const
//     {

//         assert(this->nk == p.nk);
//         assert(this->nc == p.nr);

//         Matrix result(p.nk, this->nr, p.nc);

//         for(int ik = 0; ik < this->nk; ik++)
//         {
//             result[ik] =  (*this)[ik] * p[ik];
//         }
//         return result;
//     }

//     template<typename T>
//     Matrix<T>  operator +(double s, const Matrix<T>  & p)
//     {
//         return p+s;
//     }
//     template<typename T>
//     Matrix<T>  operator *(double s, const Matrix<T>  & p)
//     {
//         return p*s;
//     }

//     template<typename T>
//     Matrix<T> operator -(double s, const Matrix<T>  & p)
//     {
//         return -p+s;
//     }
//     template<typename T>
//     bool equal_dimension(const Matrix<T>& a, const Matrix <T>& b) 
//     {
//         return a.get_nk() == b.get_nk() && a.get_nc() == b.get_nc() && a.get_nr() == b.get_nr();
//     }

    // template<typename T>
    // double  Matrix<T>::norm()
    //     {
    //         return arma::norm(data);
    //     }

//     template<typename T>
// 	Matrix<T> Matrix<T>::t() const
//     {
//         Matrix result(this->nk, this->nc, this->nr);
//         for(int ik = 0; ik < this->nk; ik++)
//         {
//             result[ik] = (*this)[ik].t();
//         }
//         return result;
//     }

//     template<typename T>
// 	Matrix<T> Matrix<T>::inv() const
//     {
//         Matrix result(*this);
//         for(int ik = 0; ik < this->nk; ik++)
//         {
//             result[ik] = arma::inv((*this)[ik]);  // Use 'result' instead of modifying '*this'
//         }
//         return result;  // Return the result, leaving the original object unchanged
//     }

//     void template<typename T>
// 	Matrix<T> Matrix<T>::print( const char* s) const
//     {
//         std::cout << s << std::endl;
//         for(int ik = 0; ik < this->nk; ik++)
//         {
//             (*this)[ik].print(std::to_string(ik));
//         }
//     }
//     void template<typename T>
// 	Matrix<T>::brief_print(const char* s) const
//     {
//         std::cout << s << std::endl;

//         for(int ik = 0; ik < this->nk; ik++)
//         {
//             (*this)[ik].brief_print(std::to_string(ik));
//         }
//     }


};
