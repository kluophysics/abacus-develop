#include "product_manifold.h"

namespace Module_Optimizer
{
    template<typename T, typename... Manifolds>
    ProductManifold<T, Manifolds...>::ProductManifold(Manifolds... manifolds) : manifolds(manifolds...) {}

    template<typename T, typename... Manifolds>
    double ProductManifold<T, Manifolds...>::metric(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldVector &xix) const
    {
        return metric_impl(x, etax, xix, std::index_sequence_for<Manifolds...>{});
    }

    template<typename T, typename... Manifolds>
    typename ProductManifold<T, Manifolds...>::ManifoldVector ProductManifold<T, Manifolds...>::projection(const ManifoldPoint &x, const ManifoldVector &etax) const
    {
        return projection_impl(x, etax, std::index_sequence_for<Manifolds...>{});
    }

    template<typename T, typename... Manifolds>
    typename ProductManifold<T, Manifolds...>::ManifoldPoint ProductManifold<T, Manifolds...>::retraction(const ManifoldPoint &x, const ManifoldVector &etax) const
    {
        return retraction_impl(x, etax, std::index_sequence_for<Manifolds...>{});
    }

    template<typename T, typename... Manifolds>
    typename ProductManifold<T, Manifolds...>::ManifoldVector ProductManifold<T, Manifolds...>::inverse_retraction(const ManifoldPoint &x, const ManifoldPoint &y) const
    {
        return inverse_retraction_impl(x, y, std::index_sequence_for<Manifolds...>{});
    }

    template<typename T, typename... Manifolds>
    typename ProductManifold<T, Manifolds...>::ManifoldVector ProductManifold<T, Manifolds...>::diff_retraction(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const
    {
        return diff_retraction_impl(x, etax, y, xix, std::index_sequence_for<Manifolds...>{});
    }

    template<typename T, typename... Manifolds>
    typename ProductManifold<T, Manifolds...>::ManifoldVector ProductManifold<T, Manifolds...>::vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const
    {
        return vector_transport_impl(x, etax, y, xix, std::index_sequence_for<Manifolds...>{});
    }

    template<typename T, typename... Manifolds>
    typename ProductManifold<T, Manifolds...>::ManifoldVector ProductManifold<T, Manifolds...>::inverse_vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xiy) const
    {
        return inverse_vector_transport_impl(x, etax, y, xiy, std::index_sequence_for<Manifolds...>{});
    }

    template<typename T, typename... Manifolds>
    int ProductManifold<T, Manifolds...>::dimension() const
    {
        return dimension_impl(std::index_sequence_for<Manifolds...>{});
    }

    template<typename T, typename... Manifolds>
    template<std::size_t... Is>
    double ProductManifold<T, Manifolds...>::metric_impl(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldVector &xix, std::index_sequence<Is...>) const
    {
        return (std::get<Is>(manifolds).metric(std::get<Is>(x), std::get<Is>(etax), std::get<Is>(xix)) + ...);
    }

    template<typename T, typename... Manifolds>
    template<std::size_t... Is>
    typename ProductManifold<T, Manifolds...>::ManifoldVector ProductManifold<T, Manifolds...>::projection_impl(const ManifoldPoint &x, const ManifoldVector &etax, std::index_sequence<Is...>) const
    {
        return std::make_tuple(std::get<Is>(manifolds).projection(std::get<Is>(x), std::get<Is>(etax))...);
    }

    template<typename T, typename... Manifolds>
    template<std::size_t... Is>
    typename ProductManifold<T, Manifolds...>::ManifoldPoint ProductManifold<T, Manifolds...>::retraction_impl(const ManifoldPoint &x, const ManifoldVector &etax, std::index_sequence<Is...>) const
    {
        return std::make_tuple(std::get<Is>(manifolds).retraction(std::get<Is>(x), std::get<Is>(etax))...);
    }

    template<typename T, typename... Manifolds>
    template<std::size_t... Is>
    typename ProductManifold<T, Manifolds...>::ManifoldVector ProductManifold<T, Manifolds...>::inverse_retraction_impl(const ManifoldPoint &x, const ManifoldPoint &y, std::index_sequence<Is...>) const
    {
        return std::make_tuple(std::get<Is>(manifolds).inverse_retraction(std::get<Is>(x), std::get<Is>(y))...);
    }

    template<typename T, typename... Manifolds>
    template<std::size_t... Is>
    typename ProductManifold<T, Manifolds...>::ManifoldVector ProductManifold<T, Manifolds...>::diff_retraction_impl(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix, std::index_sequence<Is...>) const
    {
        return std::make_tuple(std::get<Is>(manifolds).diff_retraction(std::get<Is>(x), std::get<Is>(etax), std::get<Is>(y), std::get<Is>(xix))...);
    }

    template<typename T, typename... Manifolds>
    template<std::size_t... Is>
    typename ProductManifold<T, Manifolds...>::ManifoldVector ProductManifold<T, Manifolds...>::vector_transport_impl(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix, std::index_sequence<Is...>) const
    {
        return std::make_tuple(std::get<Is>(manifolds).vector_transport(std::get<Is>(x), std::get<Is>(etax), std::get<Is>(y), std::get<Is>(xix))...);
    }

    template<typename T, typename... Manifolds>
    template<std::size_t... Is>
    typename ProductManifold<T, Manifolds...>::ManifoldVector ProductManifold<T, Manifolds...>::inverse_vector_transport_impl(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xiy, std::index_sequence<Is...>) const
    {
        return std::make_tuple(std::get<Is>(manifolds).inverse_vector_transport(std::get<Is>(x), std::get<Is>(etax), std::get<Is>(y), std::get<Is>(xiy))...);
    }

    template<typename T, typename... Manifolds>
    template<std::size_t... Is>
    int ProductManifold<T, Manifolds...>::dimension_impl(std::index_sequence<Is...>) const
    {
        return (std::get<Is>(manifolds).dimension() + ...);
    }
}

// Explicit template instantiation
template class Module_Optimizer::ProductManifold<double, Module_Optimizer::StiefelManifold<double>, 
Module_Optimizer::StiefelManifold<std::complex<double>>>;