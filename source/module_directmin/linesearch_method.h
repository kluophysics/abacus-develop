#ifndef LINESEARCH_METHOD_H
#define LINESEARCH_METHOD_H

#include <string>

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

    };

}



#endif //LINESEARCH_METHOD_H