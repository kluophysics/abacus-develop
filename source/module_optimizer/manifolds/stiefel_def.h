
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

}