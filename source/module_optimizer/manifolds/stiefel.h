#ifndef STIEFEL_MANIFOLD_H
#define STIEFEL_MANIFOLD_H

#include "manifold.h"
#include <complex>

namespace Module_Optimizer
{
    template<typename T>
    class StiefelManifold : public Manifold<T>
    {
    public:
        using typename Manifold<T>::ManifoldPoint;
        using typename Manifold<T>::ManifoldVector;

        StiefelManifold(int p, int n);

        double metric(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldVector &xix) const override;

        ManifoldVector projection(const ManifoldPoint &x, const ManifoldVector &etax) const override;

        ManifoldPoint retraction(const ManifoldPoint &x, const ManifoldVector &etax) const override;

        ManifoldVector inverse_retraction(const ManifoldPoint &x, const ManifoldPoint &y) const override;

        ManifoldVector diff_retraction(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const override;

        ManifoldVector vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const override;

        ManifoldVector inverse_vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xiy) const override;

        int dimension() const override;

    private:
        int p; // Number of rows
        int n; // Number of columns
    };
}

// #include "stiefel_manifold.cpp"

#endif // STIEFEL_MANIFOLD_H