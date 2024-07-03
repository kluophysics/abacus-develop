#ifndef ESolver_DirectMin_KS_LCAO_H
#define ESolver_DirectMin_KS_LCAO_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


#include "esolver_ks_lcao.h"

#include "module_psi/psi.h"
#include "module_hamilt_general/hamilt.h"
#include "module_directmin/directmin.h"
// #include "module_directmin/linesearch_options.h"
// #include "linesearch_options.h"

// without the occupation optimization, kluo 2024/07/02
namespace ModuleESolver
{
    // template<typename T, typename Device = psi::DEVICE_CPU>
    template<typename TK, typename TR>
    class ESolver_DirectMin_KS_LCAO: public ESolver_KS_LCAO<TK, TR>
    {
    public:
        ESolver_DirectMin_KS_LCAO()
        {
            this->classname = "ESolver_DirectMin_KS_LCAO";
            this->basisname = "LCAO";
        }
        void before_all_runners(Input & inp, UnitCell & cell) override;
        void runner(const int istep, UnitCell& cell) override;
        double cal_energy() override;
        void cal_force(ModuleBase::matrix& force) override;
        void cal_stress(ModuleBase::matrix& stress) override;
        void after_all_runners() override;

    private:

        // ModuleDirectMin::LineSearchOptions * ls_opts;
        // ModuleDirectMin::LineSearchOptions *ls_opts;
        // wavefunction coefficients
        // hamilt::Hamilt<T, Device>* p_hamilt = nullptr;

        ModuleDirectMin::DirectMin * p_directmin;


    };  // end class ESolver_DirectMin_KS_LCAO
} // namespace ModuleESolver

#endif // ESolver_DirectMin_KS_LCAO_H