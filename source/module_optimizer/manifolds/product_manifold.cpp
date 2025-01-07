#include "product_manifold.h"

namespace Module_Optimizer
{
    template<typename T, typename... Manifolds>
    ProductManifold<T, Manifolds...>::ProductManifold(Manifolds... manifolds) : manifolds(manifolds...) {}

    template<typename T, typename... Manifolds>
    double ProductManifold<T, Manifolds...>::metric(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldVector &xix) const
    {
        return std::apply([&](const auto&... ms) {
            return (ms.metric(std::get<ManifoldPoint>(x), std::get<ManifoldVector>(etax), std::get<ManifoldVector>(xix)) + ...);
        }, manifolds);
    }

    template<typename T, typename... Manifolds>
    typename ProductManifold<T, Manifolds...>::ManifoldVector ProductManifold<T, Manifolds...>::projection(const ManifoldPoint &x, const ManifoldVector &etax) const
    {
        return std::apply([&](const auto&... ms) {
            return std::make_tuple(ms.projection(std::get<ManifoldPoint>(x), std::get<ManifoldVector>(etax))...);
        }, manifolds);
    }

    template<typename T, typename... Manifolds>
    typename ProductManifold<T, Manifolds...>::ManifoldPoint ProductManifold<T, Manifolds...>::retraction(const ManifoldPoint &x, const ManifoldVector &etax) const
    {
        return std::apply([&](const auto&... ms) {
            return std::make_tuple(ms.retraction(std::get<ManifoldPoint>(x), std::get<ManifoldVector>(etax))...);
        }, manifolds);
    }

    template<typename T, typename... Manifolds>
    typename ProductManifold<T, Manifolds...>::ManifoldVector ProductManifold<T, Manifolds...>::inverse_retraction(const ManifoldPoint &x, const ManifoldPoint &y) const
    {
        return std::apply([&](const auto&... ms) {
            return std::make_tuple(ms.inverse_retraction(std::get<ManifoldPoint>(x), std::get<ManifoldPoint>(y))...);
        }, manifolds);
    }

    template<typename T, typename... Manifolds>
    typename ProductManifold<T, Manifolds...>::ManifoldVector ProductManifold<T, Manifolds...>::diff_retraction(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const
    {
        return std::apply([&](const auto&... ms) {
            return std::make_tuple(ms.diff_retraction(std::get<ManifoldPoint>(x), std::get<ManifoldVector>(etax), std::get<ManifoldPoint>(y), std::get<ManifoldVector>(xix))...);
        }, manifolds);
    }

    template<typename T, typename... Manifolds>
    typename ProductManifold<T, Manifolds...>::ManifoldVector ProductManifold<T, Manifolds...>::vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const
    {
        return std::apply([&](const auto&... ms) {
            return std::make_tuple(ms.vector_transport(std::get<ManifoldPoint>(x), std::get<ManifoldVector>(etax), std::get<ManifoldPoint>(y), std::get<ManifoldVector>(xix))...);
        }, manifolds);
    }

    template<typename T, typename... Manifolds>
    typename ProductManifold<T, Manifolds...>::ManifoldVector ProductManifold<T, Manifolds...>::inverse_vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xiy) const
    {
        return std::apply([&](const auto&... ms) {
            return std::make_tuple(ms.inverse_vector_transport(std::get<ManifoldPoint>(x), std::get<ManifoldVector>(etax), std::get<ManifoldPoint>(y), std::get<ManifoldVector>(xiy))...);
        }, manifolds);
    }

    template<typename T, typename... Manifolds>
    int ProductManifold<T, Manifolds...>::dimension() const
    {
        return std::apply([&](const auto&... ms) {
            return (ms.dimension() + ...);
        }, manifolds);
    }
}

// Explicit template instantiation
// template class 
// Module_Optimizer::ProductManifold<double, Module_Optimizer::StiefelManifold<double>, 
// Module_Optimizer::StiefelManifold<std::complex<double>>>;