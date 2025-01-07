#ifndef PRODUCT_MANIFOLD_H
#define PRODUCT_MANIFOLD_H

#include "manifold.h"
#include <tuple>
#include <utility> // For std::index_sequence and std::index_sequence_for

namespace Module_Optimizer
{
    template<typename T, typename... Manifolds>
    class ProductManifold : public Manifold<T>
    {
    public:
        using ManifoldPoint = std::tuple<typename Manifolds::ManifoldPoint...>;
        using ManifoldVector = std::tuple<typename Manifolds::ManifoldVector...>;

        ProductManifold(Manifolds... manifolds);

        double metric(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldVector &xix) const override;

        ManifoldVector projection(const ManifoldPoint &x, const ManifoldVector &etax) const override;

        ManifoldPoint retraction(const ManifoldPoint &x, const ManifoldVector &etax) const override;

        ManifoldVector inverse_retraction(const ManifoldPoint &x, const ManifoldPoint &y) const override;

        ManifoldVector diff_retraction(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const override;

        ManifoldVector vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const override;

        ManifoldVector inverse_vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xiy) const override;

        int dimension() const override;

    private:
        std::tuple<Manifolds...> manifolds;
    };
}

// #include "product_manifold.cpp"

#endif // PRODUCT_MANIFOLD_H