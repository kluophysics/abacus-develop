#include "stiefel.h"
#include <cassert>

#include "module_base/tool_quit.h"

namespace ModuleDirectMin
{

    // for Stiefel manifold

    double Stiefel::metric(const StiefelPoint & x, const StiefelVector & etax, const StiefelVector &xix)
    {

        assert(equal_dimension(x, etax) && equal_dimension(x, xix));

        arma::cx_mat term1, term2;

        double result=0.0;

        if(metric_type == CANONICAL)
        {
            for (int ik = 0; ik < x.get_nk(); ik++)
            {
                term1 = x[ik] * ((x[ik].t() * xix[ik]) / 2.0);
                term2 = xix[ik] - term1;
                arma::cx_double innerProduct = arma::trace(etax[ik].t() * term2);

                result += innerProduct.real(); // the metric in product manifold is the sum of metrics in submanifold
            }
            return result;
        }
        else if (metric_type == EUCLIDEAN)
        {
            for (int ik = 0; ik < x.get_nk(); ik++)
            {
                // result += arma::norm(A[ik].t() * B[ik]);
                arma::cx_double innerProduct = arma::trace(etax[ik].t() * xix[ik]);
                result += innerProduct.real(); // the metric in product manifold is the sum of metrics in submanifold
            } 

            return result;
        }
        else
        {
            ModuleBase::WARNING_QUIT("Stiefel::metric", "unknown metric type, please specify either EUCLIDEAN or CANONICAL!");
            return 0.0;
        }
    }

    StiefelPoint Stiefel::retraction(const StiefelPoint & x, const StiefelVector & etax)  
    {
        if(retraction_type == RT_POLAR)
        {
            // left for future change, same for RT_CAYLEY
            return projection(x, etax);
        }
        else if (retraction_type == RT_CAYLEY)
        {
            // left for future choices of vector transport
            // here only projection vector transport is done
            return projection(x, etax);
        }
        else // default to projection
        {
            // left for future choices of vector transport
            // here only projection vector transport is done
            return projection(x, etax);
        }     
    }

    StiefelVector Stiefel::vector_transport(const StiefelPoint &x, const StiefelVector &etax, const StiefelPoint &y, const StiefelVector &xix)
    {
        if(vector_transport_type == VT_PROJECTION)
        {
            return projection(x, etax);
        }

        printf("Error: VectorTransport has not been done!\n");
        return projection(x, etax);
        // return Manifold::VectorTransport(x, etax, y, xix, result);
        ;
    } 
    StiefelVector Stiefel::projection(const StiefelPoint &X, const StiefelVector &Z)
    {
        // Projection onto the tangent space to the Stiefel manifold.
        // Input: X is the foot point, matrix, and Z is the direction, matrix.
        // Output: P is the projection tangent vector, matrix.
        
        // Stiefel XZ = X.t() * Z;
        // Stiefel P = Z - 0.5 * X * (XZ + XZ.t());

        assert(equal_dimension(X, Z));

        StiefelVector XZ(Z);
        StiefelVector P(Z);

        // int nk = X.nk;
        XZ = X.t() * Z; // XZ[ik] is square matrix
        P = Z - 0.5* X *( XZ + XZ.t() );
    
        return P;
    }

    StiefelVector Stiefel::diff_retraction(const StiefelPoint &x, const StiefelVector &etax, const StiefelPoint &y, const StiefelVector &xix)
    {
        // NOTE: THIS HAS NOT BEEN DONE PROPERLY!


        // assert(equal_dimension(x, etax) && equal_dimension(y, xix)); 
        // // currently only support QR retraction 
        // Stiefel X(*this);

        // assert(X.nk == Z.nk);
        // assert(X.nr == Z.nr);
        // assert(X.nc == Z.nc);
        
        // arma::cx_mat Y, R;
        // // arma::cx_mat W;
        // arma::cx_mat U, D;

        // arma::cx_mat ZRinv, XZRinv, XXZRinv;

        // Stiefel W = X + Z;

        StiefelVector result(x);

        // for (int ik = 0; ik < X.nk; ik++)
        // {
        //     // W = X[ik] + Z[ik];
        //     arma::qr_econ(Y, R, W[ik]);
        //     // Apply the sign function to the diagonal of R
        //     R.diag() = arma::sign(arma::sign(R.diag()) + 0.5);
        //     R = arma::diagmat(R.diag()) * R;
        //     // // Modify X and R based on the sign of the diagonal Stiefels of R
        //     // Y = Y * arma::diagmat(arma::sign(arma::sign(arma::diagvec(R)) + 0.5));
        //     // result[ik] = Y;

        //     // note it uses the R matrix from retraction function
        //     ZRinv = Z [ik] * arma::inv(R);
        //     XZRinv = X[ik].t() * ZRinv;
        //     XXZRinv = X[ik] * XZRinv;
        //     U = arma::trimatl(XZRinv, -1); // Lower triangular part
        //     D = X[ik] * (U - U.t()) + ZRinv - XXZRinv;
        //     result[ik] = D;
        // }
        return result;
    }



    // for StiefelPoint



}
