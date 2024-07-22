#include "problem.h"

namespace ModuleDirectMin
{
    Problem::~Problem()
    {
        ;
    }

    Domain Problem::rgrad(const Domain& X)
    {
        int nk = X.p1.nk;

        Domain FX(X), G(X);
        FX = grad(X);

        // FX.brief_print("FX");

        Stiefel XFX(X.p1.t()*FX.p1);
        // Domain XFX( X );
        // X.t()*FX
        G.p1 = FX.p1 - X.p1 * XFX.t();
        G.p2 = X.p2;
        // for (int ik = 0; ik < nk; ++ik)
        // {
        //     XFX[ik] = X[ik].t() * FX[ik];
        //     G[ik] = FX[ik] - X[ik]*XFX[ik].t(); // See Edelman 2.53
        // }
        // Domain FX, XFX, G;

        // XFX = X.t() * FX;
        // G = FX - X*XFX.t(); // See Edelman 2.53
        return G;
    }

    Domain Problem::preconditioner(const Domain & C_in)
    {
        return C_in; // default Identity, means no preconditioner
    }

    void Problem::set_domain(Domain * domain_in)
    {
        domain = domain_in;
    }

}