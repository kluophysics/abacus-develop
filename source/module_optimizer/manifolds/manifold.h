#ifndef MANIFOLD_H
#define MANIFOLD_H

#include <armadillo>

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
}

#endif // MANIFOLD_H