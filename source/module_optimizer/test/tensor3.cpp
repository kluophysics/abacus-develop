#include "tensor3.h"

namespace Module_Optimizer
{
    // constructor
    template <typename T>
    Tensor3<T>::Tensor3()
    {
        nrows = 0;
        ncols = 0;
        nslices = 0;
        size = 0;
       arma::Cube<T>   data(nrows, ncols, nslices);
    } // constructor

    template <typename T>
    Tensor3<T>::Tensor3(unsigned int nrows, unsigned int ncols, unsigned int nslices)
    {
        this->nrows = nrows;
        this->ncols = ncols;
        this->nslices = nslices;
        size = nrows * ncols * nslices;
        arma::Cube<T>   data(nrows, ncols, nslices);
    }

    template <typename T>
    Tensor3<T>::Tensor3(const Tensor3 & p)
    {
        nrows = p.nrows;
        ncols = p.ncols;
        nslices = p.nslices;
        size = p.size;
        data = p.data;
    } // copy constructor

    template<typename T>
    bool equal_dimension(const Tensor3 <T> & a, const Tensor3<T>  & b)
    {
        return (a.nrows == b.nrows) && (a.ncols == b.ncols) && (a.nslices == b.nslices);
    }


    template<typename T>
    void Tensor3<T>::set_zeros()
    {
        data.zeros();
    }
    
    template<typename T>
    void Tensor3<T>::clear()
    {
        nrows = 0;
        ncols = 0;
        nslices = 0;
        size = 0;
        data.clear();
    }


}
