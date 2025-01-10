#include "problem.h"

namespace Module_Optimizer
{
    using ManifoldPoint = typename Manifold::ManifoldPoint;
    using ManifoldVector = typename Manifold::ManifoldVector;

    Problem::~Problem()
    {
    }

    double Problem::obj(const ManifoldPoint &x) const
    {
        // Provide a default implementation or leave it pure virtual
        return 0.0;
    }

    ManifoldVector Problem::grad(const ManifoldPoint &x) const
    {
        // Provide a default implementation or leave it pure virtual
        return ManifoldVector();
    }

    void Problem::evaluate_obj_and_grad(const ManifoldPoint &x) const
    {
        // Provide a default implementation or leave it pure virtual
    }

    ManifoldVector Problem::rie_grad(const ManifoldPoint &x) const
    {
        // Provide a default implementation or leave it pure virtual
        return ManifoldVector();
    }

// template <typename T>
// Problem<T>::~Problem()
//     {
//     }

// template <typename T>
// double Problem<T>::obj(const ManifoldPoint& x)
// {
//     // Provide a default implementation or leave it pure virtual
//     return 0.0;
// }

//     template <typename T>
//     typename Problem<T>::ManifoldVector & Problem<T>::grad(const ManifoldPoint & x) {
//         // Provide a default implementation or leave it pure virtual
//         return ManifoldVector();
//     }

//     template <typename T>
//     void Problem<T>::evaluate_obj_and_grad(const ManifoldPoint & x) {
//         // Provide a default implementation or leave it pure virtual
//     }

//     template <typename T>
//      typename Problem<T>::ManifoldVector & Problem<T>::RieGrad(const ManifoldPoint & x) {
//         // Provide a default implementation or leave it pure virtual
//         return ManifoldVector();
//     }


    // template class Problem<double>;
    // template class Problem<std::complex<double>>;
}