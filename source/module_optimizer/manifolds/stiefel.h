#ifndef STIEFEL_H
#define STIEFEL_H


#include "stiefel_def.h"
#include "manifold.h"

namespace Module_Optimizer
{

    class StiefelPoint;
    class StiefelVector;
    class Stiefel;

    class Stiefel:public Manifold
    {
        //Define the Riemannian metric: g_x(etax, xix). The default one is the Euclidean metric.
        virtual double metric(const StiefelPoint & x, const StiefelVector & etax, const StiefelVector &xix) ;

        // Compute the retraction result = R_x(etax). Default: result = x + etax
        virtual StiefelPoint retraction(const StiefelPoint & x, const StiefelVector & etax) ;
        
        // Computes the vector transport, i.e., result = \mathcal{T}_etax (xix)
        virtual StiefelVector vector_transport(const StiefelPoint &x, const StiefelVector &etax, const StiefelPoint &y, const StiefelVector &xix) ;
        



        // This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} etax
		virtual StiefelVector projection(const StiefelPoint &x, const StiefelVector &etax) ;

        // Compute the inverse retraction result = R_x^{-1} (etax). Default: result = x + etax
        virtual StiefelVector inverse_retraction(const StiefelPoint & x, const StiefelPoint & y) ;
        
        // Computes the vector transport by differentiated retraction, i.e., result = \mathcal{T}_{R_etax} (xix)
        virtual StiefelVector diff_retraction(const StiefelPoint &x, const StiefelVector &etax, const StiefelPoint &y, const StiefelVector &xix) ;


        // Computes the inverse vector transport, i.e., result = \mathcal{T}_etax^{-1} (xiy)
        virtual StiefelVector inverse_vector_transport(const StiefelPoint &x, const StiefelVector &etax,  const StiefelPoint &y, const StiefelVector &xiy) ;


        friend bool equal_dimension(const StiefelPoint & a, const StiefelPoint & b); // 
        friend bool equal_dimension(const StiefelPoint & a, const StiefelVector & b); // 
        friend bool equal_dimension(const StiefelVector & a, const StiefelPoint & b); // 
        friend bool equal_dimension(const StiefelVector & a, const StiefelVector & b); // 

        MetricType metric_type; // Riemannian metric
        RetractionType retraction_type; // retraction method
        VectorTranportType vector_transport_type; // vector transport method
    };


}


#endif //STIEFEL_H