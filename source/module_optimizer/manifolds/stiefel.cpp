#include "stiefel.h"
#include <cassert>
#include <armadillo>

#include "module_base/tool_quit.h"

namespace Module_Optimizer
{
    using ManifoldPoint = Manifold::ManifoldPoint;
    using ManifoldVector = Manifold::ManifoldVector;

    double StiefelManifold::metric(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldVector &xix) const
    {
        assert(x.n_slices == etax.n_slices && x.n_slices == xix.n_slices);
        double result = 0.0;

        for (size_t i = 0; i < x.n_slices; ++i)
        {
           arma::cx_mat xi = x.slice(i);
           arma::cx_mat etaxi = etax.slice(i);
           arma::cx_mat xixi = xix.slice(i);

            if (metric_type == CANONICAL)
            {
                result += std::real(arma::trace(etaxi.t() * (arma::mat(n, p, arma::fill::eye) - xi * xi.t() / 2) * xixi));
            }
            else if (metric_type == EUCLIDEAN)
            {
                result += std::real(arma::trace(etaxi.t() * xixi));
            }
            else
            {
                ModuleBase::WARNING_QUIT("Stiefel::metric", "unknown metric type, please specify either EUCLIDEAN or CANONICAL!");
                result = 0.0;
            }
        }
        return result;
    }

   ManifoldVector StiefelManifold::projection(const ManifoldPoint &x, const ManifoldVector &etax) const
    {
        ManifoldVector result(etax.n_rows, etax.n_cols, etax.n_slices);
        for (size_t i = 0; i < x.n_slices; ++i)
        {
            result.slice(i) = etax.slice(i) - x.slice(i) * arma::symmatu(x.slice(i).t() * etax.slice(i));
        }
        return result;
    }

    ManifoldPoint StiefelManifold::retraction(const ManifoldPoint &x, const ManifoldVector &etax) const {
        ManifoldPoint result(x.n_rows, x.n_cols, x.n_slices);
        switch (retraction_type) {
            case RT_QF: {
                for (size_t i = 0; i < x.n_slices; ++i) {
                    arma::cx_mat u, v;
                    arma::vec s;
                    arma::svd(u, s, v, x.slice(i) + etax.slice(i));
                    result.slice(i) = u * v.t();
                }
                break;
            }
            case RT_EXP: {
                for (size_t i = 0; i < x.n_slices; ++i) {
                    arma::cx_mat xi = x.slice(i);
                    arma::cx_mat etaxi = etax.slice(i);
                    arma::cx_mat expm = arma::expmat(xi * etaxi.t() - etaxi * xi.t());
                    result.slice(i) = expm * xi;
                }
                break;
            }
            case RT_CAYLEY: {
                for (size_t i = 0; i < x.n_slices; ++i) {
                    arma::cx_mat xi = x.slice(i);
                    arma::cx_mat etaxi = etax.slice(i);
                    arma::cx_mat I = arma::eye<arma::cx_mat>(p, p);
                    arma::cx_mat A = xi + etaxi;
                    arma::cx_mat B = xi - etaxi;
                    result.slice(i) = arma::solve(I + 0.5 * B, I - 0.5 * A) * xi;
                }
                break;
            }
            case RT_POLAR: {
                for (size_t i = 0; i < x.n_slices; ++i) {
                    arma::cx_mat xi = x.slice(i);
                    arma::cx_mat etaxi = etax.slice(i);
                    arma::cx_mat Y = xi + etaxi;
                    arma::cx_mat Q, R;
                    arma::qr(Q, R, Y);
                    result.slice(i) = Q;
                }
                break;
            }
            default:
                ModuleBase::WARNING_QUIT("Stiefel::retraction", "unknown retraction type!");
                break;
        }
        return result;
    }


    ManifoldVector StiefelManifold::inverse_retraction(const ManifoldPoint &x, const ManifoldPoint &y) const
    {
        ManifoldVector result(x.n_rows, x.n_cols, x.n_slices);
        for (size_t i = 0; i < x.n_slices; ++i)
        {
            result.slice(i) = y.slice(i) - x.slice(i);
        }
        return result;
    }


    ManifoldVector StiefelManifold::diff_retraction(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const
    {
        return projection(y, xix);
    }


    ManifoldVector StiefelManifold::vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const
    {
        ManifoldVector result(xix.n_rows, xix.n_cols, xix.n_slices);
         switch (vector_transport_type) {
            case VT_PROJECTION: {
                result = projection(y, xix);
                break;
            }
            case VT_PARALLELTRANSLATION: {
                result = xix;
                break;
            }
            case VT_DIFFERENTIATED: {
                result = diff_retraction(x, etax, y, xix);
                break;
            }
            case VT_CAYLEY: {
                for (size_t i = 0; i < x.n_slices; ++i) {
                    arma::cx_mat xi = x.slice(i);
                    arma::cx_mat etaxi = etax.slice(i);
                    arma::cx_mat xiy = xix.slice(i);
                    arma::cx_mat I = arma::eye<arma::cx_mat>(p, p);
                    arma::cx_mat A = xi + etaxi;
                    arma::cx_mat B = xi - etaxi;
                    result.slice(i) = arma::solve(I + 0.5 * B, I - 0.5 * A) * xiy;
                }
                break;
            }
            case VT_RIGGING: {
                for (size_t i = 0; i < x.n_slices; ++i) {
                    arma::cx_mat xi = x.slice(i);
                    arma::cx_mat etaxi = etax.slice(i);
                    arma::cx_mat xiy = xix.slice(i);
                    arma::cx_mat I = arma::eye<arma::cx_mat>(p, p);    
                    arma::cx_mat A = xi + etaxi;
                    arma::cx_mat B = xi - etaxi;
                    result.slice(i) = arma::solve(I + 0.5 * B, I - 0.5 * A) * xiy;
                }
                break;
            }
            case VT_PARALLELIZATION: {
                result = xix;
                break;
            }
            default:
                ModuleBase::WARNING_QUIT("Stiefel::vector_transport", "unknown vector transport type!");
                break;
        }
        return result;
    }
    
    ManifoldVector StiefelManifold::inverse_vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xiy) const
    {
        return diff_retraction(y, xiy, x, xiy);
    }

    
    int StiefelManifold::dimension() const
    {
        return p * n * k - (n * (n - 1) / 2) * k;
    }
}

// // Explicit template instantiation
// template class Module_Optimizer::StiefelManifold<double>;
// template class Module_Optimizer::StiefelManifold<std::complex<double>>;