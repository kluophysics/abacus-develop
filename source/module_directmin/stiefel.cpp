#include "stiefel.h"
#include "options.h"
#include <string>
#include <cassert>

namespace ModuleDirectMin
{

    // using namespace Manifold;

    Stiefel::Stiefel()
    {
        nk = 0;
        nr = 0;
        nc = 0;
        size = 0;
        psm.clear();

        metric_type = CANONICAL;
        retraction_type = QF;
        vector_transport_type = PROJECTION;

    }

    Stiefel::Stiefel(const Stiefel& var)
    {
        nk = var.nk;
        nr = var.nr;
        nc = var.nc;
        size = var.size;
        psm = var.psm;

        metric_type = CANONICAL;
        retraction_type = QF;
        vector_transport_type = PROJECTION;
    }

    Stiefel::Stiefel(int k, int r, int c)
    {
        nk = k;
        nr = r;
        nc = c;
        size = k*r*c;
        psm.resize(k);
        for (int ik = 0; ik < nk; ik++)
        {
            // psm[ik] = ModuleBase::ComplexMatrix(nr, nc);
            psm[ik] = arma::cx_mat(nr, nc);
        }
    }


    void Stiefel::switch_metric(MetricType met)
    {
        metric_type = met;
    }
    void Stiefel::switch_retraction(RetractionType retr)
    {
        retraction_type = retr;
    }
    void Stiefel::switch_vector_transport(VectorTranportType vectran)
    {
        vector_transport_type = vectran;
    }


    void Stiefel::resize(int k, int r, int c)
    {
        psm.resize(k);

        size = k*r*c;
        nk = k;
        nr = r;
        nc = c;
        for(int ik = 0; ik < k; ik++)
        {
            // psm[ik] = ModuleBase::ComplexMatrix(r, c);
            psm[ik] = arma::cx_mat(nr, nc, arma::fill::zeros);
        }
    }

    void Stiefel::reset()
    {
        for( int ik = 0; ik < nk; ik++)
        {
            psm[ik].zeros();
        }
    }

    void Stiefel::clear()
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

    // Stiefel& Stiefel::operator-()
    // {
    //     for(int ik = 0; ik < nk; ik++)
    //     {
    //         *this[ik] = -*this[ik];
    //     }
    //     return *this;
    // }

    // Stiefel& Stiefel::operator+()
    // {
    //     // for(int ik = 0; ik < p.nk; ik++)
    //     // {
    //     //     *this[ik] = *this[ik];
    //     // }
    //     return *this;
    // }

    arma::cx_mat Stiefel::operator[](const int k) const
    {
        assert( k >= 0 && k < nk);
        return psm[k];
    }

    arma::cx_mat & Stiefel::operator[](const int k) 
    {
        assert( k >= 0 && k < nk);
        return psm[k];
    }


    Stiefel Stiefel::operator-() const
    {
        // Stiefel result(p);

        for(int ik = 0; ik < this->nk; ik++)
        {
            (*this)[ik] = -(*this)[ik] ;
        }
        return *this;

    }


    Stiefel& Stiefel::operator=(const Stiefel&p)
    {
        this->resize(p.nk, p.nr, p.nc);

        for(int ik = 0; ik < p.nk; ik++)
        {
            (*this)[ik] = p[ik];
        }
        return *this;
    }

    Stiefel& Stiefel::operator-=(const Stiefel&p)
    {
        this->resize(p.nk, p.nr, p.nc);

        for(int ik = 0; ik < p.nk; ik++)
        {
            (*this)[ik] -= p[ik];
        }
        return *this;
    }

    Stiefel& Stiefel::operator+=(const Stiefel&p)
    {
        this->resize(p.nk, p.nr, p.nc);

        for(int ik = 0; ik < p.nk; ik++)
        {
            (*this)[ik] += p[ik];
        }
        return *this;
    }

    Stiefel& Stiefel::operator*=(const Stiefel&p)
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

    Stiefel Stiefel::operator +(double s) const
    {
        for(int ik = 0; ik < this->nk; ik++)
        {
            (*this)[ik] += s;
        }
        return *this;
    }

    Stiefel Stiefel::operator *(double s) const
    {
        for(int ik = 0; ik < this->nk; ik++)
        {
            (*this)[ik] *= s;
        }
        return *this;
    }

    Stiefel Stiefel::operator -(double s) const
    {
        for(int ik = 0; ik < this->nk; ik++)
        {
            (*this)[ik] -= s;
        }
        return *this;
    }

    Stiefel Stiefel::operator+( const Stiefel & p) const
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
    Stiefel Stiefel::operator-( const Stiefel & p) const
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

    Stiefel Stiefel::operator*( const Stiefel & p) const
    {

        assert(this->nk == p.nk);
        assert(this->nc == p.nr);

        Stiefel result(p.nk, this->nr, p.nc);

        for(int ik = 0; ik < this->nk; ik++)
        {
            result[ik] =  (*this)[ik] * p[ik];
        }
        return result;
    }


    Stiefel operator +(double s, const Stiefel & p)
    {
        return p+s;
    }
    Stiefel operator *(double s, const Stiefel & p)
    {
        return p*s;
    }
    Stiefel operator -(double s, const Stiefel & p)
    {
        return -p+s;
    }

    double Stiefel::norm()
    {
        double result = 0.0;
        for(int ik = 0; ik < this->nk; ik++)
        {
            result += arma::norm( (*this)[ik] );
        }
        return result/this->nk;
    }

    Stiefel Stiefel::t() const
    {
        // Stiefel result(*this);
        Stiefel result(this->nk, this->nc, this->nr);
        for(int ik = 0; ik < this->nk; ik++)
        {
            result[ik] = (*this)[ik].t();
            // result[ik] = (*this)[ik].t();
        }
        return result;
    }

    Stiefel Stiefel::inv() const
    {
        Stiefel result(*this);
        for(int ik = 0; ik < this->nk; ik++)
        {
            result[ik] = arma::inv((*this)[ik]);  // Use 'result' instead of modifying '*this'
        }
        return result;  // Return the result, leaving the original object unchanged
    }

    void Stiefel::print( const char* s) const
    {
        std::cout << s << std::endl;
        for(int ik = 0; ik < this->nk; ik++)
        {
            (*this)[ik].print(std::to_string(ik));
        }
    }
    void Stiefel::brief_print(const char* s) const
    {
        std::cout << s << std::endl;

        for(int ik = 0; ik < this->nk; ik++)
        {
            (*this)[ik].brief_print(std::to_string(ik));
        }
    }




    // const Stiefel operator+(const Stiefel&p,const   Stiefel& q)
    // {
    //     assert(q.nk == p.nk);
    //     assert(q.nr == p.nr);
    //     assert(q.nc == p.nc);

    //     Stiefel result(p);

    //     for(int ik = 0; ik < p.nk; ik++)
    //     {
    //         result[ik] = q[ik] + p[ik];
    //     }
    //     return result;
    // }

    // const Stiefel operator-(const Stiefel&p,const   Stiefel& q)
    // {
    //     assert(q.nk == p.nk);
    //     assert(q.nr == p.nr);
    //     assert(q.nc == p.nc);

    //     Stiefel result(p);

    //     for(int ik = 0; ik < p.nk; ik++)
    //     {
    //         result[ik] =  p[ik] -q[ik] ;
    //     }
    //     return result;
    // }

    // Stiefel operator*(const Stiefel&p,const   Stiefel& q)
    // {
    //     assert(q.nk == p.nk);
    //     assert(p.nc == q.nr);

    //     Stiefel result(p.nk, p.nr, q.nc);

    //     for(int ik = 0; ik < p.nk; ik++)
    //     {
    //         result[ik] = p[ik] * q[ik];
    //     } 
    //     return result;
    // }



    // Stiefel operator/(const Stiefel&p, const   Stiefel& q)
    // {
    //     assert(q.nk == p.nk);
    //     assert(q.nr == p.nr);
    //     assert(q.nc == p.nc);
    //     Stiefel result(p.nk, p.nr, q.nc);

    //     for(int ik = 0; ik < p.nk; ik++)
    //     {
    //         result[ik] = p[ik] / q[ik];
    //     }
    //     return result;
    // }

    // Stiefel operator+(const Stiefel&p,  double s)
    // {
    //     Stiefel result;

    //     for(int ik = 0; ik < p.nk; ik++)
    //     {
    //         result[ik] = s+p[ik] ;
    //     }
    //     return result;
    // }
    // Stiefel operator+(double s, const Stiefel&p )
    // {
    //     Stiefel result(p+s);

    //     // for(int ik = 0; ik < p.nk; ik++)
    //     // {
    //     //     result[ik] = s+p[ik] ;
    //     // }
    //     return result;
    // }

    // Stiefel operator-(const Stiefel&p,  double s)
    // {
    //     Stiefel result;

    //     for(int ik = 0; ik < p.nk; ik++)
    //     {
    //         result[ik] = p[ik] - s ;
    //     }
    //     return result;
    // }

    // Stiefel operator-( double s, const Stiefel&p )
    // {
    //     Stiefel result(s - p);

    //     // for(int ik = 0; ik < p.nk; ik++)
    //     // {
    //     //     result[ik] = s - p[ik] ;
    //     // }
    //     return result;
    // }
    // Stiefel operator*(const Stiefel&p,  double s)
    // {
    //     Stiefel result;

    //     for(int ik = 0; ik < p.nk; ik++)
    //     {
    //         result[ik] = s*p[ik] ;
    //     }
    //     return result;
    // }

    // Stiefel operator*( double s, const Stiefel&p)
    // {

    //     Stiefel result(p);

    //     for(int ik = 0; ik < p.nk; ik++)
    //     {
    //         result[ik] = s*p[ik] ;
    //     }
    //     return result;
    // }


    // Stiefel proj(const Stiefel& X, const Stiefel& Z)
    // {
    //     // Projection onto the tangent space to the Stiefel manifold.
    //     // Input: X is the foot point, matrix, and Z is the direction, matrix.
    //     // Output: P is the projection tangent vector, matrix.
        
    //     // Stiefel XZ = X.t() * Z;
    //     // Stiefel P = Z - 0.5 * X * (XZ + XZ.t());
        
    //     assert(X.size == Z.size);
    //     int nk = X.nk;

    //     Stiefel XZ(X);
    //     Stiefel P(X);

    //     // for (int ik = 0; ik < nk; ++ik)
    //     // {
    //     //     XZ[ik] = X[ik].t() * Z[ik]; // XZ[ik] is square matrix
    //     //     P[ik] = Z[ik] - 0.5* X[ik]*( XZ[ik] + XZ[ik].t() );
    //     // }

    //     XZ = X.t() * Z; // XZ[ik] is square matrix
    //     P = Z - 0.5* X *( XZ + XZ.t() );
        

    //     return P;
    // }

    // Stiefel vectran(const Stiefel& X, const Stiefel& Z)
    // {
    //     // if (vectran_choice == 'qr')
    //     // left for future choices of vector transport
    //     // here only projection vector transport is done
    //     return proj(X, Z);
    // }


    // Stiefel retraction(const Stiefel& X, const Stiefel& Z)
    // {
    //     assert(X.nk == Z.nk);
    //     assert(X.nr == Z.nr);
    //     assert(X.nc == Z.nc);

    //     // arma::cx_mat W;
    //     Stiefel W = X + Z;
    //     Stiefel result(X);

    //     for (int ik = 0; ik < X.nk; ik++)
    //     {
    //         arma::cx_mat Y, R;
    //         arma::qr_econ(Y, R, W[ik]);
    //         // Apply the sign function to the diagonal of R
    //         R.diag() = arma::sign(arma::sign(R.diag()) + 0.5);
    //         // // Modify X and R based on the sign of the diagonal Stiefels of R
    //         Y = Y * arma::diagmat(arma::sign(arma::sign(arma::diagvec(R)) + 0.5));

    //         result[ik] = Y;
    //     }
    //     return result;
    // }

    // Stiefel diff_retraction(const Stiefel& X, const Stiefel& Z)
    // {
    //     assert(X.nk == Z.nk);
    //     assert(X.nr == Z.nr);
    //     assert(X.nc == Z.nc);
        
    //     arma::cx_mat Y, R;
    //     // arma::cx_mat W;
    //     arma::cx_mat U, D;

    //     arma::cx_mat ZRinv, XZRinv, XXZRinv;

    //     Stiefel W = X + Z;

    //     Stiefel result(X);

    //     for (int ik = 0; ik < X.nk; ik++)
    //     {
    //         // W = X[ik] + Z[ik];
    //         arma::qr_econ(Y, R, W[ik]);
    //         // Apply the sign function to the diagonal of R
    //         R.diag() = arma::sign(arma::sign(R.diag()) + 0.5);
    //         R = arma::diagmat(R.diag()) * R;
    //         // // Modify X and R based on the sign of the diagonal Stiefels of R
    //         // Y = Y * arma::diagmat(arma::sign(arma::sign(arma::diagvec(R)) + 0.5));
    //         // result[ik] = Y;

    //         // note it uses the R matrix from retraction function
    //         ZRinv = Z [ik] * arma::inv(R);
    //         XZRinv = X[ik].t() * ZRinv;
    //         XXZRinv = X[ik] * XZRinv;
    //         U = arma::trimatl(XZRinv, -1); // Lower triangular part
    //         D = X[ik] * (U - U.t()) + ZRinv - XXZRinv;
    //         result[ik] = D;
    //     }
    //     return result;
    // }

    // double metric(const Stiefel& X, const Stiefel& A, const Stiefel& B)
    // {
    //     assert(A.nk == B.nk);
    //     assert(A.nr == B.nr);
    //     assert(A.nc == B.nc);

    //     assert(X.nk == A.nk);
    //     assert(X.nr == A.nr);
    //     assert(X.nc == A.nc);

    //     arma::cx_mat term1, term2;

    //     // (X.t()*A).brief_print("Y'*A");
    //     // (X.t()*B).brief_print("Y'*B");

    //     double result=0.0;

    //     bool use_canonical_metric = true;
    //     if(use_canonical_metric)

    //         for (int ik = 0; ik < X.nk; ik++)
    //         {
    //             term1 = X[ik] * ((X[ik].t() * B[ik]) / 2.0);
    //             term2 = B[ik] - term1;
    //             arma::cx_double innerProduct = arma::trace(A[ik].t() * term2);

    //             result += innerProduct.real(); // the metric in product manifold is the sum of metrics in submanifold
    //         }
    //     // std::cout << "Canonical metric=" << result << std::endl;

       
    //     if(false)
    //     {
    //          result = 0.0;
    //         for (int ik = 0; ik < X.nk; ik++)
    //         {
    //             // result += arma::norm(A[ik].t() * B[ik]);
    //             arma::cx_double innerProduct = arma::trace(A[ik].t() * B[ik]);
    //             result += innerProduct.real(); // the metric in product manifold is the sum of metrics in submanifold
    //         }
    //     std::cout << "Euclidean metric=" << result << std::endl;
    //     }        

    //     return result;
    // }

    // bool is_orthogonal(const Stiefel& X)
    // {
    //     bool result = false;

    //     for(int ik = 0; ik < X.nk; ik++)
    //     {
            
    //     }
    //     return false;
    // }

    // Stiefel vectran(const Stiefel& X, const Stiefel& Z)
    // {
    //     return proj(X, Z);
    //     // return Stiefel();
    // }

    // Stiefel diff_retraction(const Stiefel& X, const Stiefel& Z)
    // {
    //     return Stiefel();
    // }

    // double metric(const Stiefel& X, const Stiefel& A, Stiefel& B)
    // {
    //     return 0.0;
    // }


    Stiefel Stiefel::projection(const Stiefel& Z)
    {
        // Projection onto the tangent space to the Stiefel manifold.
        // Input: X is the foot point, matrix, and Z is the direction, matrix.
        // Output: P is the projection tangent vector, matrix.
        
        // Stiefel XZ = X.t() * Z;
        // Stiefel P = Z - 0.5 * X * (XZ + XZ.t());

        Stiefel X(*this);
        Stiefel XZ(X);
        Stiefel P(X);

        assert(X.size == Z.size);
        // int nk = X.nk;
        XZ = X.t() * Z; // XZ[ik] is square matrix
        P = Z - 0.5* X *( XZ + XZ.t() );
    
        return P;

    }

    Stiefel Stiefel::vector_transport(const Stiefel& Z)
    {
        if(vector_transport_type == DIFFERENTIATED)
        {
            return this->projection(Z);
        }
        else if (vector_transport_type == CAYLEYVT)
        {
            // left for future choices of vector transport
            // here only projection vector transport is done
            return this->projection(Z);
        }
        else
        {
            // left for future choices of vector transport
            // here only projection vector transport is done
            return this->projection(Z);
        }     
        // if (vectran_choice == 'qr')
        // left for future choices of vector transport
        // here only projection vector transport is done
        // return proj(X, Z);
    }


    Stiefel Stiefel::retraction(const Stiefel& Z)
    {
        // currently only support QR retraction 
        Stiefel X(*this);

        assert(X.nk == Z.nk);
        assert(X.nr == Z.nr);
        assert(X.nc == Z.nc);

        // arma::cx_mat W;
        Stiefel W = X + Z;
        Stiefel result(X);

        if(retraction_type == EXP)
        {
            return result; 

        }
        else if(retraction_type == CAYLEY)
        {
            return result;
        }
        else if(retraction_type == POLAR)
        {
            return result;
        }
        else // STIE_QF is default
        {
            for (int ik = 0; ik < X.nk; ik++)
            {
                arma::cx_mat Y, R;
                arma::qr_econ(Y, R, W[ik]);
                // Apply the sign function to the diagonal of R
                R.diag() = arma::sign(arma::sign(R.diag()) + 0.5);
                // // Modify X and R based on the sign of the diagonal Stiefels of R
                Y = Y * arma::diagmat(arma::sign(arma::sign(arma::diagvec(R)) + 0.5));

                result[ik] = Y;
            }
            return result;

        }

    }

    // Stiefel Stiefel::diff_retraction(const Stiefel& Z)
    // {        
    //     // currently only support QR retraction 
    //     Stiefel X(*this);

    //     assert(X.nk == Z.nk);
    //     assert(X.nr == Z.nr);
    //     assert(X.nc == Z.nc);
        
    //     arma::cx_mat Y, R;
    //     // arma::cx_mat W;
    //     arma::cx_mat U, D;

    //     arma::cx_mat ZRinv, XZRinv, XXZRinv;

    //     Stiefel W = X + Z;

    //     Stiefel result(X);

    //     for (int ik = 0; ik < X.nk; ik++)
    //     {
    //         // W = X[ik] + Z[ik];
    //         arma::qr_econ(Y, R, W[ik]);
    //         // Apply the sign function to the diagonal of R
    //         R.diag() = arma::sign(arma::sign(R.diag()) + 0.5);
    //         R = arma::diagmat(R.diag()) * R;
    //         // // Modify X and R based on the sign of the diagonal Stiefels of R
    //         // Y = Y * arma::diagmat(arma::sign(arma::sign(arma::diagvec(R)) + 0.5));
    //         // result[ik] = Y;

    //         // note it uses the R matrix from retraction function
    //         ZRinv = Z [ik] * arma::inv(R);
    //         XZRinv = X[ik].t() * ZRinv;
    //         XXZRinv = X[ik] * XZRinv;
    //         U = arma::trimatl(XZRinv, -1); // Lower triangular part
    //         D = X[ik] * (U - U.t()) + ZRinv - XXZRinv;
    //         result[ik] = D;
    //     }
    //     return result;
    // }
    double Stiefel::metric(const Stiefel& A, const Stiefel& B)
    {
        if(metric_type == CANONICAL)
        {
            return canonical_metric(A,B);
        }
        else if(metric_type == EUCLIDEAN)
        {
            return euclidean_metric(A,B);
        }
        else
        {
            return 0.0;
        }
    }

    double Stiefel::canonical_metric(const Stiefel& A, const Stiefel& B)
    {
        Stiefel X(*this);

        assert(A.nk == B.nk);
        assert(A.nr == B.nr);
        assert(A.nc == B.nc);

        assert(X.nk == A.nk);
        assert(X.nr == A.nr);
        assert(X.nc == A.nc);

        arma::cx_mat term1, term2;

        // (X.t()*A).brief_print("Y'*A");
        // (X.t()*B).brief_print("Y'*B");

        double result=0.0;

        for (int ik = 0; ik < X.nk; ik++)
        {
            term1 = X[ik] * ((X[ik].t() * B[ik]) / 2.0);
            term2 = B[ik] - term1;
            arma::cx_double innerProduct = arma::trace(A[ik].t() * term2);

            result += innerProduct.real(); // the metric in product manifold is the sum of metrics in submanifold
        }
        return result;
       
  

    }

    double Stiefel::euclidean_metric(const Stiefel& A, const Stiefel& B)
    {
        Stiefel X(*this);

        assert(A.nk == B.nk);
        assert(A.nr == B.nr);
        assert(A.nc == B.nc);

        assert(X.nk == A.nk);
        assert(X.nr == A.nr);
        assert(X.nc == A.nc);

        arma::cx_mat term1, term2;

        double result=0.0;
       
        for (int ik = 0; ik < X.nk; ik++)
        {
            // result += arma::norm(A[ik].t() * B[ik]);
            arma::cx_double innerProduct = arma::trace(A[ik].t() * B[ik]);
            result += innerProduct.real(); // the metric in product manifold is the sum of metrics in submanifold
        } 

        return result;
    }
    bool Stiefel::is_orthogonal()
    {
        Stiefel X(*this);

        bool result = false;

        for(int ik = 0; ik < X.nk; ik++)
        {
            
        }
        return false;
    }
};
