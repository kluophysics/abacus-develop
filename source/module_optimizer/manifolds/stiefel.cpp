#include "stiefel.h"

namespace Module_Optimizer
{
    template<typename T>
    StiefelManifold<T>::StiefelManifold(int p, int n) : p(p), n(n) {}

    template<typename T>
    double StiefelManifold<T>::metric(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldVector &xix) const
    {
        // if constexpr (std::is_same<T, std::complex<double>>::value)
        // if (std::is_same<T, std::complex<double>>::value)
        // {
        //     return std::real(arma::cdot(arma::vectorise(etax), arma::vectorise(xix)));
        // }
        // else
        // {
        //     return std::real(arma::dot(arma::vectorise(etax), arma::vectorise(xix)));
        // }
        return std::real(arma::dot(arma::vectorise(etax), arma::vectorise(xix)));

    }

    template<typename T>
    typename StiefelManifold<T>::ManifoldVector StiefelManifold<T>::projection(const ManifoldPoint &x, const ManifoldVector &etax) const
    {
        return etax - x * arma::symmatu(x.t() * etax);
    }

    template<typename T>
    typename StiefelManifold<T>::ManifoldPoint StiefelManifold<T>::retraction(const ManifoldPoint &x, const ManifoldVector &etax) const
    {
        arma::Mat<T> u, v;
        arma::Col<typename arma::get_pod_type<T>::result> s;
        arma::svd(u, s, v, x + etax);
        return u * v.t();
    }

    template<typename T>
    typename StiefelManifold<T>::ManifoldVector StiefelManifold<T>::inverse_retraction(const ManifoldPoint &x, const ManifoldPoint &y) const
    {
        return y - x;
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
        return p * n - (n * (n - 1)) / 2;
    }
}

// Explicit template instantiation
template class Module_Optimizer::StiefelManifold<double>;
template class Module_Optimizer::StiefelManifold<std::complex<double>>;