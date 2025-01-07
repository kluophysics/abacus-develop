#ifndef MANIFOLD_H
#define MANIFOLD_H

#include <armadillo>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Module_Optimizer
{
    template<typename T>
    class Manifold
    {
    public:
        using ManifoldPoint = arma::Mat<T>;
        using ManifoldVector = arma::Mat<T>;

        virtual ~Manifold() = default;

        virtual double metric(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldVector &xix) const = 0;
        virtual ManifoldVector projection(const ManifoldPoint &x, const ManifoldVector &etax) const = 0;
        virtual ManifoldPoint retraction(const ManifoldPoint &x, const ManifoldVector &etax) const = 0;
        virtual ManifoldVector inverse_retraction(const ManifoldPoint &x, const ManifoldPoint &y) const = 0;
        virtual ManifoldVector diff_retraction(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const = 0;
        virtual ManifoldVector vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const = 0;
        virtual ManifoldVector inverse_vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xiy) const = 0;
        virtual int dimension() const = 0;
    };

    template<typename... Manifolds>
    class ProductManifold : public Manifold<typename std::common_type<typename Manifolds::ManifoldPoint...>::type::elem_type>
    {
    public:
        using T = typename std::common_type<typename Manifolds::ManifoldPoint...>::type::elem_type;
        using ManifoldPoint = std::tuple<typename Manifolds::ManifoldPoint...>;
        using ManifoldVector = std::tuple<typename Manifolds::ManifoldVector...>;

        ProductManifold(Manifolds... manifolds) : manifolds_(manifolds...) {}

        double metric(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldVector &xix) const override
        {
            return metric_impl(x, etax, xix, std::index_sequence_for<Manifolds...>{});
        }

        ManifoldVector projection(const ManifoldPoint &x, const ManifoldVector &etax) const override
        {
            return projection_impl(x, etax, std::index_sequence_for<Manifolds...>{});
        }

        ManifoldPoint retraction(const ManifoldPoint &x, const ManifoldVector &etax) const override
        {
            return retraction_impl(x, etax, std::index_sequence_for<Manifolds...>{});
        }

        ManifoldVector inverse_retraction(const ManifoldPoint &x, const ManifoldPoint &y) const override
        {
            return inverse_retraction_impl(x, y, std::index_sequence_for<Manifolds...>{});
        }

        ManifoldVector diff_retraction(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const override
        {
            return diff_retraction_impl(x, etax, y, xix, std::index_sequence_for<Manifolds...>{});
        }

        ManifoldVector vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const override
        {
            return vector_transport_impl(x, etax, y, xix, std::index_sequence_for<Manifolds...>{});
        }

        ManifoldVector inverse_vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xiy) const override
        {
            return inverse_vector_transport_impl(x, etax, y, xiy, std::index_sequence_for<Manifolds...>{});
        }

        int dimension() const override
        {
            return dimension_impl(std::index_sequence_for<Manifolds...>{});
        }

    private:
        std::tuple<Manifolds...> manifolds_;

        template<std::size_t... Is>
        double metric_impl(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldVector &xix, std::index_sequence<Is...>) const
        {
            return (std::get<Is>(manifolds_).metric(std::get<Is>(x), std::get<Is>(etax), std::get<Is>(xix)) + ...);
        }

        template<std::size_t... Is>
        ManifoldVector projection_impl(const ManifoldPoint &x, const ManifoldVector &etax, std::index_sequence<Is...>) const
        {
            return std::make_tuple(std::get<Is>(manifolds_).projection(std::get<Is>(x), std::get<Is>(etax))...);
        }

        template<std::size_t... Is>
        ManifoldPoint retraction_impl(const ManifoldPoint &x, const ManifoldVector &etax, std::index_sequence<Is...>) const
        {
            return std::make_tuple(std::get<Is>(manifolds_).retraction(std::get<Is>(x), std::get<Is>(etax))...);
        }

        template<std::size_t... Is>
        ManifoldVector inverse_retraction_impl(const ManifoldPoint &x, const ManifoldPoint &y, std::index_sequence<Is...>) const
        {
            return std::make_tuple(std::get<Is>(manifolds_).inverse_retraction(std::get<Is>(x), std::get<Is>(y))...);
        }

        template<std::size_t... Is>
        ManifoldVector diff_retraction_impl(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix, std::index_sequence<Is...>) const
        {
            return std::make_tuple(std::get<Is>(manifolds_).diff_retraction(std::get<Is>(x), std::get<Is>(etax), std::get<Is>(y), std::get<Is>(xix))...);
        }

        template<std::size_t... Is>
        ManifoldVector vector_transport_impl(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix, std::index_sequence<Is...>) const
        {
            return std::make_tuple(std::get<Is>(manifolds_).vector_transport(std::get<Is>(x), std::get<Is>(etax), std::get<Is>(y), std::get<Is>(xix))...);
        }

        template<std::size_t... Is>
        ManifoldVector inverse_vector_transport_impl(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xiy, std::index_sequence<Is...>) const
        {
            return std::make_tuple(std::get<Is>(manifolds_).inverse_vector_transport(std::get<Is>(x), std::get<Is>(etax), std::get<Is>(y), std::get<Is>(xiy))...);
        }

        template<std::size_t... Is>
        int dimension_impl(std::index_sequence<Is...>) const
        {
            return (std::get<Is>(manifolds_).dimension() + ...);
        }
    };
}

#endif // MANIFOLD_H