#include "problem.h"

namespace ModuleDirectMin
{
    Problem::~Problem()
    {
        ;
    }

    Composite Problem::rgrad(const Composite& X)
    {
        int nk = X.p1.nk;

        Composite FX(X), G(X);
        FX = grad(X);

        // FX.brief_print("FX");

        Stiefel XFX(X.p1.t()*FX.p1);
        // Composite XFX( X );
        // X.t()*FX
        G.p1 = FX.p1 - X.p1 * XFX.t();
        G.p2 = X.p2;
        // for (int ik = 0; ik < nk; ++ik)
        // {
        //     XFX[ik] = X[ik].t() * FX[ik];
        //     G[ik] = FX[ik] - X[ik]*XFX[ik].t(); // See Edelman 2.53
        // }
        // Composite FX, XFX, G;

        // XFX = X.t() * FX;
        // G = FX - X*XFX.t(); // See Edelman 2.53
        return G;
    }

    Composite Problem::preconditioner(const Composite & C_in)
    {
        return C_in; // default Identity, means no preconditioner
    }

    void Problem::set_domain(Composite * domain_in)
    {
        domain = domain_in;
    }

}