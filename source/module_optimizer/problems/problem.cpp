#include "problem.h"

namespace Module_Optimizer
{
    template <typename T>
    Problem<T>::~Problem()
    {
    };



    template class Problem<double>;
    template class Problem<std::complex<double>>;
}