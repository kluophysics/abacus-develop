#ifndef FLEXMAT_H
#define FLEXMAT_H

#include <armadillo>

namespace ModuleDirectMin
{
    // This FlexMat is implemented with the help of armadillo library https://arma.sourceforge.net.
    // In the future, it can be replaced by a more efficient implementation
    // or a local variant that could be parallalized 
    // Note: this new version directly used Cube format for the product of manifolds
    //       in comparison to the original vector of matrices type.

    // since the matrix for a particular k is often referenced, it seems that 
    // let k be the slice, then the matrix uses a continious segment of memory, which should be fast
    // (correnct me if I'm wrong!!!)
    class FlexMat
    {
        public:

        int nk; // number of manifolds, each manifold in this case is of same size
        int nr; // number of rows for each manifold
        int nc; // number of columns for each manifold

        int size; // nk*nr*nc

        arma::cx_cube data; // actual data

        
    };

}

#endif // FLEXMAT_H