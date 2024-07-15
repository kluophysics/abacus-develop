#ifndef STIEFEL_H
#define STIEFEL_H

#define StiefelVector StiefelPoint
#include "matrix_vector.h"

#include "manifold.h"

namespace ModuleDirectMin
{
    /* Riemannian Metric for the Stiefel manifold:
	Eucldean: g_x(etax, xix) = \trace(etax^T xix);
	Canonical: g_x(etax, xix) = \trace(etax^T (I_n - x x^T / 2) xix); */
	enum MetricType { 
        EUCLIDEAN, 
        CANONICAL, 
        MetricTypeLength
    };

	/*Retraction for the Stiefel manifold
	QF: qf retraction defined in [AMS2008, (4.8)]
	POLAR: polar based retraction defined in [AMS2008, (4.7)]
	EXP: The exponential mapping
	CAYLEYR: the Cayley transform in [Zhu2016]
	[AMS2008]P.-A. Absil, R. Mahony, and R. Sepulchre. Optimization algorithms on matrix manifolds.
	Princeton University Press, Princeton, NJ, 2008.
	[HGA2015]:Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
	SIAM Journal on Optimization, 25(3):1660?685,2015.
	[Hua2013]:W. Huang. Optimization algorithms on Riemannian manifolds with applications.
	PhD thesis, Florida State University, Department of Mathematics, 2013.
	[Zhu2016]: Xiaojing Zhu, A Riemannian conjugate gradient method for optimization on the Stiefel Manifold */
    enum RetractionType {
        QF,
        EXP,
        CAYLEY,
        POLAR,
        RetractionTypeLength
    };
    // VectorTransport should  use the corresponding retraction method 

	/*Vector transport for the Stiefel manifold
	PARALLELIZATION: Vector transport by parallelization, See [HAG2015, Section 2.3.1]
	RIGGING: Vector transport by rigging, See [HAG2015, Section 2.3.2]
	PARALLELTRANSLATION: parallel translation
	CAYLEYVT: the vector transport based on Cayley transform. [Zhu2016]
	[HAG2015]:W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trust-region method.
	Mathematical Programming, 150(2):179?16, February 2015
	[Zhu2016]: Xiaojing Zhu, A Riemannian conjugate gradient method for optimization on the Stiefel Manifold */
    enum VectorTranportType {
        PROJECTION,
        PARALLELTRANSLATION,
        DIFFERENTIATED,
        CAYLEYVT,
        RIGGING,
        PARALLELIZATION,
        VectorTranportTypeLength
    };

    class StiefelPoint;
    class StiefelVector;

    class Stiefel: public Manifold 
    {
    public:

        //Define the Riemannian metric: g_x(etax, xix). The default one is the Euclidean metric.
        virtual double metric(const StiefelPoint & x, const StiefelVector & etax, const StiefelVector &xix) const;

        // This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} etax
		virtual StiefelVector & projection(const StiefelPoint &x, const StiefelVector &etax) const;

        // Compute the retraction result = R_x(etax). Default: result = x + etax
        virtual StiefelPoint & retraction(const StiefelPoint & x, const StiefelVector & etax) const;
        
        // Compute the inverse retraction result = R_x^{-1} (etax). Default: result = x + etax
        virtual StiefelVector & inverse_retraction(const StiefelPoint & x, const StiefelPoint & y) const;
        
        // Computes the vector transport by differentiated retraction, i.e., result = \mathcal{T}_{R_etax} (xix)
        virtual StiefelVector &diff_retraction(const StiefelPoint &x, const StiefelVector &etax, const StiefelPoint &y, const StiefelVector &xix) const;

        // Computes the vector transport, i.e., result = \mathcal{T}_etax (xix)
        virtual StiefelVector &vector_transport(const StiefelPoint &x, const StiefelVector &etax, const StiefelVector &xix) const;
        
        // Computes the inverse vector transport, i.e., result = \mathcal{T}_etax^{-1} (xiy)
        virtual StiefelVector &inverse_vector_transport(const StiefelPoint &x, const StiefelVector &etax, const StiefelVector &xiy) const;

        inline void switch_metric(MetricType metric_type_in) 
        {
            metric_type = metric_type_in;
        };
        inline void switch_retraction(RetractionType retraction_type_in)
        {
            retraction_type = retraction_type_in;
        };
        inline void switch_vector_transport(VectorTranportType vector_transport_type_in)
        {
            vector_transport_type = vector_transport_type_in;
        };

        MetricType metric_type; // Riemannian metric
        RetractionType retraction_type; // retraction method
        VectorTranportType vector_transport_type; // vector transport method
    };

    // StiefelPoint inherits from MatrixVector class and so does the StiefelVector
    class StiefelPoint: public ManifoldPoint, public MatrixVector
    {
        bool is_orthogonal(); // check for orthogonal
    };

    // // StiefelVector inherits from MatrixVector class 
    // class StiefelVector: public ManifoldVector , public MatrixVector
    // {
    //     bool is_orthogonal(); // check for orthogonal
    // };
}

#endif //STIEFEL_H