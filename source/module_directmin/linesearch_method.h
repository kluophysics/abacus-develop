#ifndef LINESEARCH_METHOD_H
#define LINESEARCH_METHOD_H

#include <string>
#include <vector>

namespace ModuleDirectMin
{

    class LineSearchMethodBase
    {
    public:
        LineSearchMethodBase();
        virtual ~LineSearchMethodBase()=0;
        ;
        std::string method_name; // name of method: conjugate gradient
        bool verbose; // verbosity level


        // f1: function value at x1, f2: function value at x2
        double f1, f2;

        
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
    };

}



#endif //LINESEARCH_METHOD_H