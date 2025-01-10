#ifndef OptimizerLSBase_H
#define OptimizerLSBase_H

#include <list>

#include "optimizer_base.h"
#include "optimizer_ls_options.h"

namespace Module_Optimizer
{
    // LS is short for line search, used as a parent class for line search algorithms, 
    // such as steepest descent, conjugate gradient, quasi-Newton, etc.
    class OptimizerLSBase : public OptimizerBase
    {

        
    public:
        // the main function to optimize the problem
        virtual void optimize();

        // set default params for line search
        virtual void set_default_params();

        // update params for line search using LSOptions *opt_in
        virtual void update_params(LSOptions *opt_in);
        
        // a pure function to be overidden by derived classes
        virtual void get_search_direction()=0;
        
        // do a linesearch with the step size according to the condition type
        virtual void do_line_search();
        
        // update one iteration, some algorithms need to update some information. For example,
		// quasi-Newton methods need to update the Hessian approximation and nonlinear conjugate gradient
		// needs to update the search direction. They are done in the following function
        virtual void update_data()=0;

        // compute the initial step size as a guess
        virtual void init_step_guess();


		//Evaluate the cost function
        //      phi(stepsize) = f(R_{x_1}(stepsize * d1))
        virtual double phi();

        //Evaluate the derivative of cost function h, 
        //       h'(stepsize) = \frac{d}{d stepsize} f(R_{x_1}(stepsize * d1))
        virtual double dphi();


        //The strong Wolfe condition [NW06 Algorithm 3.5]
		// [NW06] : J.Nocedal and S.J.Wright.Numerical optimization.Springer, second edition, 2006
        virtual void StrongWolfe(); 

		/*The Armijo-Goldstein condition.[DS83 Algorithm A6.3.1] combined with nonmontone line search
		[DS83] : J.E.Dennis and R.B.Schnabel.Numerical methods for unconstrained optimization and nonlinear equations.Springer, New Jersey, 1983*/
		virtual void Armijo();


		/*The weak Wolfe condition [DS83 Algorithm A6.3.1mod]
		[DS83] : J.E.Dennis and R.B.Schnabel.Numerical methods for unconstrained optimization and nonlinear equations.Springer, New Jersey, 1983*/
		virtual void Wolfe();


        // print information 
        virtual void print_info();


        double step_size; // The stepsize
        double step_size_old; // The old steps
        double initial_length; // inital stepsize at an iteration

        double f_previous; // function value at previous step
        

        double initial_slope_pre; // intial slope
        double initial_slope; // current slope
        double new_slope; // new slope
		std::list<double> pre_funs; // Store a few computed function values for nonmonotonic line search*/



        // algorithm-related variables:
        // x1: current iterate, x2: next iterate
        ManifoldPoint x1, x2;
        //gf1: gradient at x1,  gf2: gradient at x2
        ManifoldVector gf1, gf2; 
        // d1 is the search direction. d2 is stepsize * d1.
        ManifoldVector d1, d2;
        // f1: function value at x1, f2: function value at x2
        double f1, f2;

		// ManifoldPoint Pgf1, Pgf2;	/*Pgf1: preconditioned gradient at x1, Pgf2: preconditioned gradient at x2, used in RCG and LRBFGS, LRTRSR1*/

        // ngf0: the norm of gradient at the initial iterate
        // ngf1: the norm of the gradient at x1, 
        // ngf2: the norm of the gradient at x2
        double ngf0, ngf1, ngf2; 

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


    public:
        int numPreFuns;

        ConditionType condition_type;
        double LS_alpha;;
        double LS_beta;
        double LS_c1;
        double LS_c2; 
        double initial_step_size;
        double max_step_size;
        double min_step_size;
        double final_step_size;
        double gtol;
        double ftol;
        double tolerance;
        double accuracy;


    protected:


    private:
        void zoom(double x1, double fx1, double slope1, double x2, double fx2);

    };

}

#endif // OptimizerLSBase_H