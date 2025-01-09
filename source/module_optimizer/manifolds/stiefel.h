#ifndef STIEFEL_MANIFOLD_H
#define STIEFEL_MANIFOLD_H

#include "manifold.h"
#include <complex>
#include <armadillo>

namespace Module_Optimizer
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
	RT_QF: qf retraction defined in [AMS2008, (4.8)]
	RT_POLAR: polar based retraction defined in [AMS2008, (4.7)]
	RT_EXP: The exponential mapping
	CAYLEYR: the Cayley transform in [Zhu2016]
	[AMS2008]P.-A. Absil, R. Mahony, and R. Sepulchre. Optimization algorithms on matrix manifolds.
	Princeton University Press, Princeton, NJ, 2008.
	[HGA2015]:Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
	SIAM Journal on Optimization, 25(3):1660?685,2015.
	[Hua2013]:W. Huang. Optimization algorithms on Riemannian manifolds with applications.
	PhD thesis, Florida State University, Department of Mathematics, 2013.
	[Zhu2016]: Xiaojing Zhu, A Riemannian conjugate gradient method for optimization on the Stiefel Manifold */
    enum RetractionType {
        RT_QF,
        RT_EXP,
        RT_CAYLEY,
        RT_POLAR,
        RetractionTypeLength
    };
    // VectorTransport should  use the corresponding retraction method 

	/*Vector transport for the Stiefel manifold
	VT_PARALLELIZATION: Vector transport by parallelization, See [HAG2015, Section 2.3.1]
	VT_RIGGING: Vector transport by rigging, See [HAG2015, Section 2.3.2]
	VT_PARALLELTRANSLATION: parallel translation
	VT_CAYLEY: the vector transport based on Cayley transform. [Zhu2016]
	[HAG2015]:W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trust-region method.
	Mathematical Programming, 150(2):179?16, February 2015
	[Zhu2016]: Xiaojing Zhu, A Riemannian conjugate gradient method for optimization on the Stiefel Manifold */
    enum VectorTranportType {
        VT_PROJECTION,
        VT_PARALLELTRANSLATION,
        VT_DIFFERENTIATED,
        VT_CAYLEY,
        VT_RIGGING,
        VT_PARALLELIZATION,
        VectorTranportTypeLength
    };

    
    class StiefelManifold : public Manifold
    {
    public:
        // using typename Manifold<T>::ManifoldPoint;
        // using typename Manifold<T>::ManifoldVector;

        // StiefelManifold(int p, int n);
        StiefelManifold(int nr, int nc, int num_manifolds) : 
            p(nr), n(nc), k(num_manifolds), metric_type(EUCLIDEAN), retraction_type(RT_QF), vector_transport_type(VT_PROJECTION) {}

        double metric(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldVector &xix) const override;

        ManifoldVector projection(const ManifoldPoint &x, const ManifoldVector &etax) const override;

        ManifoldPoint retraction(const ManifoldPoint &x, const ManifoldVector &etax) const override;

        ManifoldVector inverse_retraction(const ManifoldPoint &x, const ManifoldPoint &y) const override;

        ManifoldVector diff_retraction(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const override;

        ManifoldVector vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xix) const override;

        ManifoldVector inverse_vector_transport(const ManifoldPoint &x, const ManifoldVector &etax, const ManifoldPoint &y, const ManifoldVector &xiy) const override;

        int dimension() const override;
        inline void switch_metric(MetricType metric_type_in) { metric_type = metric_type_in; }
        inline void switch_retraction(RetractionType retraction_type_in) { retraction_type = retraction_type_in; }
        inline void switch_vector_transport(VectorTranportType vector_transport_type_in) { vector_transport_type = vector_transport_type_in; }
        
    private:
        int p; // Number of rows
        int n; // Number of columns
        int k; // Number of Stiefel manifolds in the product
        MetricType metric_type; // Riemannian metric
        RetractionType retraction_type; // retraction method
        VectorTranportType vector_transport_type; // vector transport method
    };
}

// #include "stiefel_manifold.cpp"

#endif // STIEFEL_MANIFOLD_H