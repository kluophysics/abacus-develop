#ifndef OPTIM_TYPES_H
#define OPTIM_TYPES_H


namespace Module_Optim
{
	/*Specify what information will be output in the algorithm.
	The value should be assigned to the member variable: "VerbosityLevel" */
    enum VerbosityLevel { LOW, MEDIUM, HIGH, COMPLETE, VerbosityLevelLength};

    /*The algorithm is stopped when a value (specified by ther parameter) is less than the "Tolerance" (a member variable)
	The value should be assigned to the member variable: "StopCriterion" and the applicable values are
	FUNC_REL: |f_k - f_{k+1}| / max(|f_k|, 1)
    FUNC_ABS: |f_k - f_{k+1}| 
	GRAD_ABS: \|gf_k\|
	GRAD_REL: \|gf_k\| / \|gf_0\|*/
    enum StopCriterion { FUNC_REL, FUNC_ABS, GRAD_ABS, GRAD_REL,  StopCriterionLength};

    /*Provide two types of Optimization types
    LS:  line search algorithm
    TR: trust-region algorithm
    */
    enum OptimizerType { LS, TR};

}


#endif // OPTIM_TYPES_H