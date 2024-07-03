#include "esolver_directmin_ks_lcao.h"
// #include "module_directmin/linesearch_options.cpp"

namespace ModuleESolver
{
    template<typename TK, typename TR>
    void ESolver_DirectMin_KS_LCAO<TK, TR>::before_all_runners(Input &inp, UnitCell &cell)
    {
        // Do all the stuff (mostly necessary) that ESolver_KS_LCAO<TK, TR> 
        // does first.
        ESolver_KS_LCAO<TK, TR>::before_all_runners(inp, cell);


        // then initialize all needed parameters that are needed for 
        // the direct minmization part
        // ls_opts = new ModuleDirectMin::LineSearchOptions(inp);

        // ls_opts()



        
        // ModuleBase::TITLE("ESolver_DirectMin_KS_LCAO", "Init");
        // GlobalV::ofs_running << "Inside ESolver_DirectMin_KS_LCAO";
        // // dm->Initialize();
        // dm = new ModuleDirectMin::DirectMin(inp, cell);
        // // std::cout << "--------------------------------"
		// // << "after " << "ESolver_DirectMin_KS_LCAO<TK, TR>::Init()" << std::endl;
        // dm->initialize(inp, cell);
    }

    template<typename TK, typename TR>
    // template<typename T, typename Device>
    // void ESolver_DirectMin_KS_LCAO<T, Device>::Run(const int istep, UnitCell &cell)
    void ESolver_DirectMin_KS_LCAO<TK, TR>::runner(const int istep, UnitCell &cell)
    {
        // std::cout << "--------------------------------"
		// << "before " << "ESolver_DirectMin_KS_LCAO<TK, TR>::Run()" << std::endl;

        // dm -> run(istep, cell);
        // std::cout << "--------------------------------"
		// << "after " << "ESolver_DirectMin_KS_LCAO<TK, TR>::Run()" << std::endl;

    }

    template<typename TK, typename TR>
    double ESolver_DirectMin_KS_LCAO<TK, TR>::cal_energy()
    {
        return 0.0;
    }

    template<typename TK, typename TR>
    void ESolver_DirectMin_KS_LCAO<TK, TR>::cal_force(ModuleBase::matrix& force)
    {
        ;
    }

    template<typename TK, typename TR>
    void ESolver_DirectMin_KS_LCAO<TK, TR>::cal_stress(ModuleBase::matrix& stress)
    {
        ;
    }

    template<typename TK, typename TR>
    void ESolver_DirectMin_KS_LCAO<TK, TR>::after_all_runners()
    {
        ;
    }

    // don't forget to include templates lastly...
    template class ESolver_DirectMin_KS_LCAO<double, double>;
    template class ESolver_DirectMin_KS_LCAO<std::complex<double>, double>;
    template class ESolver_DirectMin_KS_LCAO<std::complex<double>, std::complex<double>>;

} // namespace ModuleESolver