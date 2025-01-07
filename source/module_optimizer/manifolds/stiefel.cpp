#include "stiefel.h"
#include <cassert>
#include <armadillo>

#include "module_base/tool_quit.h"

namespace Module_Optimizer
{
    // template<typename T>
    // StiefelManifold<T>::StiefelManifold(int p, int n) : p(p), n(n) {}

    template<typename T>
    double StiefelManifold<T>::metric(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldVector &xix) const
    {
        assert(x.n_slices == etax.n_slices && x.n_slices == xix.n_slices);
        double result = 0.0;

        for (size_t i = 0; i < x.n_slices; ++i)
        {
            arma::Mat<T> xi = x.slice(i);
            arma::Mat<T> etaxi = etax.slice(i);
            arma::Mat<T> xixi = xix.slice(i);

            if (metric_type == CANONICAL)
            {
                result += std::real(arma::trace(etaxi.t() * (arma::Mat<T>(n, p, arma::fill::eye) - xi * xi.t() / 2) * xixi));
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

    template<typename T>
    typename StiefelManifold<T>::ManifoldVector StiefelManifold<T>::projection(const ManifoldPoint &x, const ManifoldVector &etax) const
    {
        ManifoldVector result(etax.n_rows, etax.n_cols, etax.n_slices);
        for (size_t i = 0; i < x.n_slices; ++i)
        {
            result.slice(i) = etax.slice(i) - x.slice(i) * arma::symmatu(x.slice(i).t() * etax.slice(i));
        }
        return result;
    }

    template<typename T>
    typename StiefelManifold<T>::ManifoldPoint StiefelManifold<T>::retraction(const ManifoldPoint &x, const ManifoldVector &etax) const
    {
        ManifoldPoint result(x.n_rows, x.n_cols, x.n_slices);
        for (size_t i = 0; i < x.n_slices; ++i)
        {
            arma::Mat<T> u, v;
            arma::Col<typename arma::get_pod_type<T>::result> s;
            arma::svd(u, s, v, x.slice(i) + etax.slice(i));
            result.slice(i) = u * v.t();
        }
        return result;
    }

    template<typename T>
    typename StiefelManifold<T>::ManifoldVector StiefelManifold<T>::inverse_retraction(const ManifoldPoint &x, const ManifoldPoint &y) const
    {
        ManifoldVector result(x.n_rows, x.n_cols, x.n_slices);
        for (size_t i = 0; i < x.n_slices; ++i)
        {
            result.slice(i) = y.slice(i) - x.slice(i);
        }
        return result;
    }

    template<typename T>
    typename StiefelManifold<T>::ManifoldVector StiefelManifold<T>::diff_retraction(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const
    {
        return projection(y, xix);
    }

    template<typename T>
    typename StiefelManifold<T>::ManifoldVector StiefelManifold<T>::vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const
    {
        return diff_retraction(x, etax, y, xix);
    }

    template<typename T>
    typename StiefelManifold<T>::ManifoldVector StiefelManifold<T>::inverse_vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xiy) const
    {
        return diff_retraction(y, xiy, x, xiy);
    }

    template<typename T>
    int StiefelManifold<T>::dimension() const
    {
        return p * n * k - (n * (n - 1) / 2) * k;
    }
}

// Explicit template instantiation
template class Module_Optimizer::StiefelManifold<double>;
template class Module_Optimizer::StiefelManifold<std::complex<double>>;