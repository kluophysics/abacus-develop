
#ifndef MANIFOLD_H
#define MANIFOLD_H


namespace Module_Optimizer
{
    class ManifoldPoint;
    class ManifoldVector;
    class Manifold;


    class ManifoldPoint 
    {

    };

    class ManifoldVector
    {
        ;
    };

    class Manifold 
    {
        //Define the Riemannian metric: g_x(etax, xix). The default one is the Euclidean metric.
        virtual double metric(const ManifoldPoint & x, const ManifoldVector & etax, const ManifoldVector &xix) ;

        // This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} etax
		virtual ManifoldVector  projection(const ManifoldPoint &x, const ManifoldVector &etax) ;

        // Compute the retraction result = R_x(etax). Default: result = x + etax
        virtual ManifoldPoint  retraction(const ManifoldPoint & x, const ManifoldVector & etax) ;
        
        // Compute the inverse retraction result = R_x^{-1} (etax). Default: result = x + etax
        virtual ManifoldVector  inverse_retraction(const ManifoldPoint & x, const ManifoldPoint & y) ;
        
        // Computes the vector transport by differentiated retraction, i.e., result = \mathcal{T}_{R_etax} (xix)
        virtual ManifoldVector diff_retraction(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) ;

        // Computes the vector transport, i.e., result = \mathcal{T}_etax (xix)
        virtual ManifoldVector vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) ;
        
        // Computes the inverse vector transport, i.e., result = \mathcal{T}_etax^{-1} (xiy)
        virtual ManifoldVector inverse_vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xiy) ;
    };

}

#endif //MANIFOLD_H