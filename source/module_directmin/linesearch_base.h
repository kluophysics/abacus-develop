#ifndef LINESEARCH_BASE_H
#define LINESEARCH_BASE_H

#include <string>
#include <vector>
#include <list>

#include "problem.h"
#include "composite.h"
#include "linesearch_options.h"

namespace ModuleDirectMin
{
    /* Linesearch status. It is an output argument and users don't need to assign this enumerate to any member variable.
	NOCURVATURE: the second Wolfe condition is not satisfied
	MINSTEPSIZE: line search algorithm reaches the minimum stepsize
	MAXSTEPSIZE: line search algorithm reaches the maximum stepsize
	NONEXACT: exact line search algorithm does not find a point satisfying the inner stopping criterion
	SUCCESS: line search algorithm succeeds in finding a point satisfying the line search condition
	*/
    enum LineSearchStatus{ 
        NOCURVATURE, 
        MINSTEPSIZE, 
        MAXSTEPSIZE, 
        // NONEXACT, 
        LSERROR, 
        SUCCESS, 
        LineSearchStatusLength };

	/*Initial step size in line search algorithm.
	ONESTEP: t0 = one 
	BBSTEP: t0 = g(s, s) / g(s, y), s is the difference of consecutive iterates and y is the difference of the
			gradients at consecutie iterates.
	QUADINT: t0 = [(3.60), NW06]
	QUADINTMOD: t0 = [page 60, NW06]
	[NW06]: J. Nocedal and S. J. Wright. Numerical optimization. Springer, second edition, 2006
	*/
    enum InitStepsizeType { 
        ONESTEP, 
        BBSTEP, 
        QUADINT, 
        QUADINTMOD, 
        EXTRBBSTEP, 
        InitStepsizeTypeLength };

    enum ConditionType {ARMIJO,
                        WOLFE,
                        STRONG_WOLFE,
                        ConditionTypeLength
                       };

    class LineSearchBase
    {
    public:
        LineSearchBase();
        virtual ~LineSearchBase()=0;
        ;

    protected:
        // set default params for line search
        virtual void set_default_parameters();
        
        // update params for line search using LineSearchOptions *ls_opt_in
        virtual void update_parameters(LineSearchOptions * ls_opt_in);

        // a pure function to be overidden by derived classes
        virtual void get_search_direction() = 0;

        // do a linesearch with the step size according to the condition type
        virtual void do_line_search();

        // iterate a step
        virtual void iterate();

		//Evaluate the cost function
        //      phi(stepsize) = f(R_{x_1}(stepsize * \eta))
        // virtual double phi();

        //Evaluate the derivative of cost function h, 
        //       dphi(stepsize) = \frac{d}{d stepsize} f(R_{x_1}(stepsize * \eta))
        // virtual double dphi();

        // Evaluate the cost function and the derivative at the same time.
        // phi(stepsize) = f(R_{x_1}(stepsize * \eta))
        // dphi(stepsize) = \frac{d}{d stepsize} f(R_{x_1}(stepsize * \eta))
        virtual void evaluate_phi_and_dphi(double *phi, double *dphi);



        //The strong Wolfe condition [NW06 Algorithm 3.5]
		// [NW06] : J.Nocedal and S.J.Wright.Numerical optimization.Springer, second edition, 2006
        virtual void StrongWolfe(); 

		/*The Armijo-Goldstein condition.[DS83 Algorithm A6.3.1] combined with nonmontone line search
		[DS83] : J.E.Dennis and R.B.Schnabel.Numerical methods for unconstrained optimization and nonlinear equations.Springer, New Jersey, 1983*/
		virtual void Armijo();


		/*The weak Wolfe condition [DS83 Algorithm A6.3.1mod]
		[DS83] : J.E.Dennis and R.B.Schnabel.Numerical methods for unconstrained optimization and nonlinear equations.Springer, New Jersey, 1983*/
		virtual void Wolfe();
        

        Problem * prob;

        // parameters for this class

        std::string method_name; // name of method: conjugate gradient
        bool verbose; // verbosity level
        LineSearchOptions * ls_options;
        ConditionType condition_type; // The condition type
        
        LineSearchStatus LS_status;

        // algorithm-related variables:
        // x1: current iterate, x2: next iterate
        Composite x1, x2;
        // gf1: gradient of current iterate, gf2: gradient of next iterate
        Composite gf1, gf2;
        /*Pgf1: preconditioned gradient at x1, Pgf2: preconditioned gradient at x2, used in RCG and RBFGS*/
        Composite Pgf1, Pgf2;

        /*In Line search-based methods, eta1 is the search direction. eta2 is stepsize * eta1.*/
        Composite eta1, eta2; 


        // f1: function value at x1, f2: function value at x2
        double f1, f2;

        double step_size; // The stepsize


        double step_size_old; // The old steps
        double initial_length; // inital stepsize at an iteration

        double f_previous; // function value at previous step
        

        double initial_slope_pre; // intial slope
        double initial_slope; // current slope
        double new_slope; // new slope
		std::list<double> pre_funs; // Store a few computed function values for nonmonotonic line search*/



        // for debug

		int iter; /*iteration number*/

        int nf; // number of function evaluations
        int ng; // number of gradient evaluations
        int nR; // number of retraction evaluations
        int nV; // number of vector transport
        int nH; // number of action of Hessian

        std::vector<double> time_series; // computational time series
        std::vector<double> fun_series;  // function series
        std::vector<double> grad_series; // gradient series
        std::vector<double> step_series; // stepsize series


    private:
        /*The function used in the strong Wolfe condition. See p.62 Algorithm 3.6 in [NW06]
		[NW06] : J.Nocedal and S.J.Wright.Numerical optimization.Springer, second edition, 2006 */
        void zoom(double x1, double fx1, double slope1, double x2, double fx2, double slope2);


        
    };

}



#endif //LINESEARCH_BASE_H