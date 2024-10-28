#include "domain.h"
#include <cassert>

namespace ModuleDirectMin
{


    // DomainPoint Domain::retraction(const DomainPoint& dp,
    //                                const DomainVector& dv)                               
    // {
    //     DomainPoint result(dp);
    //     result.x = Stiefel::retraction(dp.x, dv.x);
    //     return result;
    // }
    // DomainVector Domain::diff_retraction(const DomainPoint& dp1,
    //                                      const DomainVector dv1,
    //                                      const DomainPoint& dp2,
    //                                      const DomainVector& dv2) 
    // {
    //     DomainVector result(dp2);
    //     result.x = Stiefel::diff_retraction(dp1.x, dv1.x, dp2.x, dv2.x);
    //     return result;
    // }
    double Domain::norm(const DomainPoint& dp,
                          DomainVector& dv1,
                          const DomainVector& dv2) 
    {
        double total_metric = 0.0;
        assert( equal_dimension(dp, dv1));
        assert( equal_dimension(dv1,dv2));

        // metric contribution from the occupation vector
        if(explicit_occupation_flag)
        {
            for(int ik=0; ik < dp.get_nk(); ++ik)
            {
                total_metric += arma::norm(dp.occ.occ_vector[ik]);
            }
        }
        // metric contribution from the Stiefel variable
        total_metric += Stiefel::metric(dp.x, dv1.x, dv2.x);
        return total_metric;
    }
    DomainPoint::DomainPoint(const DomainPoint& dp) 
    {
        x = dp.x;
        occ = dp.occ;
    }

    DomainPoint::DomainPoint(StiefelPoint& x_in, const Occupation& occ_in) 
    {
        x = x_in;
        occ = occ_in;
    }

    void Domain::apply_async_scaling(DomainVector * dv, double beta) 
    {
        (*dv).occ = beta * (*dv).occ;
    }

    bool equal_dimension(const DomainPoint& a, const DomainPoint& b) 
    {
        return a.occ.nk == b.occ.nk && a.occ.norb == b.occ.norb && equal_dimension(a.x, b.x);
    }
    // bool equal_dimension(const DomainPoint& a, const DomainVector & b) 
    // {
    //     return a.occ.nk == b.occ.nk && a.occ.norb == b.occ.norb && equal_dimension(a.x, b.x);
    // }
    //     bool equal_dimension(const DomainVector& a, const DomainVector& b) 
    // {
    //     return a.occ.nk == b.occ.nk && a.occ.norb == b.occ.norb && equal_dimension(a.x, b.x);
    // }
    } // namespace ModuleDirectMin