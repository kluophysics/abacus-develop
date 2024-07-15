#include "manifold.h"


namespace ModuleDirectMin
{
    double Manifold::metric(const ManifoldPoint & x, const ManifoldVector & etax, const ManifoldVector &xix) const
    {
        ;
    }

    ManifoldVector & Manifold::projection(const ManifoldPoint &x, const ManifoldVector &etax) const
    {
        ;
    }

    ManifoldPoint & Manifold::retraction(const ManifoldPoint & x, const ManifoldVector & etax) const
    {
        ;
    }
    ManifoldVector & Manifold::inverse_retraction(const ManifoldPoint & x, const ManifoldPoint & y) const
    {
        ;
    }
    ManifoldVector & Manifold::diff_retraction(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const
    {
        ;
    }
    ManifoldVector & Manifold::vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldVector &xix) const
    {
        ;
    }
    ManifoldVector & Manifold::inverse_vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldVector &xiy) const
    {
        ;
    }

}