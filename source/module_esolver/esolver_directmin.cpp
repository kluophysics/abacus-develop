#include "esolver_directmin.h"

namespace ModuleESolver
{
    template<typename TK, typename TR>
    // template<typename T, typename Device>
    // void ESolver_DirectMin<T, Device>::Init(Input &inp, UnitCell &cell)
    void ESolver_DirectMin<TK, TR>::before_all_runners(Input &inp, UnitCell &cell)
    {
        // ESolver_FP::Init(inp, cell);
        
        // ModuleBase::TITLE("ESolver_DirectMin", "Init");
        // GlobalV::ofs_running << "Inside ESolver_DirectMin";
        // // dm->Initialize();
        // dm = new ModuleDirectMin::DirectMin(inp, cell);
        // // std::cout << "--------------------------------"
		// // << "after " << "ESolver_DirectMin<TK, TR>::Init()" << std::endl;
        // dm->initialize(inp, cell);
    }

    template<typename TK, typename TR>
    // template<typename T, typename Device>
    // void ESolver_DirectMin<T, Device>::Run(const int istep, UnitCell &cell)
    void ESolver_DirectMin<TK, TR>::runner(const int istep, UnitCell &cell)
    {
        // std::cout << "--------------------------------"
		// << "before " << "ESolver_DirectMin<TK, TR>::Run()" << std::endl;

        // dm -> run(istep, cell);
        // std::cout << "--------------------------------"
		// << "after " << "ESolver_DirectMin<TK, TR>::Run()" << std::endl;

    }

    template<typename TK, typename TR>
    double ESolver_DirectMin<TK, TR>::cal_energy()
    {
        return 0.0;
    }

    template<typename TK, typename TR>
    void ESolver_DirectMin<TK, TR>::cal_force(ModuleBase::matrix& force)
    {
        ;
    }

    template<typename TK, typename TR>
    void ESolver_DirectMin<TK, TR>::cal_stress(ModuleBase::matrix& stress)
    {
        ;
    }

    template<typename TK, typename TR>
    void ESolver_DirectMin<TK, TR>::after_all_runners()
    {
        ;
    }

    // don't forget to include templates lastly...
    template class ESolver_DirectMin<double, double>;
    template class ESolver_DirectMin<std::complex<double>, double>;
    template class ESolver_DirectMin<std::complex<double>, std::complex<double>>;

} // namespace ModuleESolver