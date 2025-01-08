#include "problem.h"

namespace Module_Optimizer
{

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


    template class Problem<double>;
    template class Problem<std::complex<double>>;
}