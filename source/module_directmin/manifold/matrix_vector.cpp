#include "matrix_vector.h"
#include <string>
#include <cassert>

namespace ModuleDirectMin
{

    MatrixVector::MatrixVector()
    {
        nk = 1;
        nr = 1;
        nc = 1;
        size = 1;
        psm.resize(nk);
    }

    MatrixVector::MatrixVector(const MatrixVector& var)
    {
        nk = var.nk;
        nr = var.nr;
        nc = var.nc;
        size = var.size;
        psm = var.psm;
    }

    MatrixVector::MatrixVector(int k, int r, int c)
    {
        nk = k;
        nr = r;
        nc = c;
        size = k*r*c;
        psm.resize(k);
        for (int ik = 0; ik < nk; ik++)
        {
            psm[ik] = arma::cx_mat(nr, nc);
        }
    }

    void MatrixVector::resize(int k, int r, int c)
    {
        psm.resize(k);

        size = k*r*c;
        nk = k;
        nr = r;
        nc = c;
        for(int ik = 0; ik < k; ik++)
        {
            psm[ik] = arma::cx_mat(nr, nc, arma::fill::zeros);
        }
    }

    void MatrixVector::reset()
    {
        for( int ik = 0; ik < nk; ik++)
        {
            psm[ik].zeros();
        }
    }

    void MatrixVector::clear()
    {
        nk = 0;
        nr = 0;
        nc = 0;
        size = 0;
        for( int ik = 0; ik < nk; ik++)
        {
            psm[ik].clear();
        }
    }


    arma::cx_mat MatrixVector::operator[](const int k) const
    {
        assert( k >= 0 && k < nk);
        return psm[k];
    }

    arma::cx_mat & MatrixVector::operator[](const int k) 
    {
        assert( k >= 0 && k < nk);
        return psm[k];
    }


    MatrixVector MatrixVector::operator-() const
    {
        // MatrixVector result(p);

        for(int ik = 0; ik < this->nk; ik++)
        {
            (*this)[ik] = -(*this)[ik] ;
        }
        return *this;

    }


    MatrixVector& MatrixVector::operator=(const MatrixVector&p)
    {
        this->resize(p.nk, p.nr, p.nc);

        for(int ik = 0; ik < p.nk; ik++)
        {
            (*this)[ik] = p[ik];
        }
        return *this;
    }

    MatrixVector& MatrixVector::operator-=(const MatrixVector&p)
    {
        this->resize(p.nk, p.nr, p.nc);

        for(int ik = 0; ik < p.nk; ik++)
        {
            (*this)[ik] -= p[ik];
        }
        return *this;
    }

    MatrixVector& MatrixVector::operator+=(const MatrixVector&p)
    {
        this->resize(p.nk, p.nr, p.nc);

        for(int ik = 0; ik < p.nk; ik++)
        {
            (*this)[ik] += p[ik];
        }
        return *this;
    }

    MatrixVector& MatrixVector::operator*=(const MatrixVector&p)
    {
        assert(this->nk == p.nk);
        assert(this->nc == p.nr);

        this->resize(p.nk, this->nr, p.nc);

        for(int ik = 0; ik < p.nk; ik++)
        {
            (*this)[ik] *= p[ik];
        }
        return *this;
    }

    MatrixVector MatrixVector::operator +(double s) const
    {
        for(int ik = 0; ik < this->nk; ik++)
        {
            (*this)[ik] += s;
        }
        return *this;
    }

    MatrixVector MatrixVector::operator *(double s) const
    {
        for(int ik = 0; ik < this->nk; ik++)
        {
            (*this)[ik] *= s;
        }
        return *this;
    }

    MatrixVector MatrixVector::operator -(double s) const
    {
        for(int ik = 0; ik < this->nk; ik++)
        {
            (*this)[ik] -= s;
        }
        return *this;
    }

    MatrixVector MatrixVector::operator+( const MatrixVector & p) const
    {
        assert(this->nk == p.nk);
        assert(this->nc == p.nc);
        assert(this->nr == p.nr);

        for(int ik = 0; ik < this->nk; ik++)
        {
            (*this)[ik] += p[ik];
        }
        return *this;
    }
    MatrixVector MatrixVector::operator-( const MatrixVector & p) const
    {
        assert(this->nk == p.nk);
        assert(this->nc == p.nc);
        assert(this->nr == p.nr);

        for(int ik = 0; ik < this->nk; ik++)
        {
            (*this)[ik] -= p[ik];
        }
        return *this;
    }

    MatrixVector MatrixVector::operator*( const MatrixVector & p) const
    {

        assert(this->nk == p.nk);
        assert(this->nc == p.nr);

        MatrixVector result(p.nk, this->nr, p.nc);

        for(int ik = 0; ik < this->nk; ik++)
        {
            result[ik] =  (*this)[ik] * p[ik];
        }
        return result;
    }


    MatrixVector operator +(double s, const MatrixVector & p)
    {
        return p+s;
    }
    MatrixVector operator *(double s, const MatrixVector & p)
    {
        return p*s;
    }
    MatrixVector operator -(double s, const MatrixVector & p)
    {
        return -p+s;
    }

    double MatrixVector::norm()
    {
        double result = 0.0;
        for(int ik = 0; ik < this->nk; ik++)
        {
            result += arma::norm( (*this)[ik] );
        }
        return result/this->nk;
    }

    MatrixVector MatrixVector::t() const
    {
        MatrixVector result(this->nk, this->nc, this->nr);
        for(int ik = 0; ik < this->nk; ik++)
        {
            result[ik] = (*this)[ik].t();
        }
        return result;
    }

    MatrixVector MatrixVector::inv() const
    {
        MatrixVector result(*this);
        for(int ik = 0; ik < this->nk; ik++)
        {
            result[ik] = arma::inv((*this)[ik]);  // Use 'result' instead of modifying '*this'
        }
        return result;  // Return the result, leaving the original object unchanged
    }

    void MatrixVector::print( const char* s) const
    {
        std::cout << s << std::endl;
        for(int ik = 0; ik < this->nk; ik++)
        {
            (*this)[ik].print(std::to_string(ik));
        }
    }
    void MatrixVector::brief_print(const char* s) const
    {
        std::cout << s << std::endl;

        for(int ik = 0; ik < this->nk; ik++)
        {
            (*this)[ik].brief_print(std::to_string(ik));
        }
    }

};
