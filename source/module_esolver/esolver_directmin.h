#ifndef ESOLVER_DIRECTMIN_H
#define ESOLVER_DIRECTMIN_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


#include "esolver_fp.h"

#include "module_psi/psi.h"
#include "module_hamilt_general/hamilt.h"
#include "module_directmin/directmin.h"

namespace ModuleESolver
{
    // template<typename T, typename Device = psi::DEVICE_CPU>
    template<typename TK, typename TR>
    class ESolver_DirectMin: public ESolver_FP
    {
    public:
        ESolver_DirectMin()
        {
            classname = "ESolver_DirectMin";
        }
        virtual void before_all_runners(Input & inp, UnitCell & cell) override;
        virtual void runner(const int istep, UnitCell& cell) override;
        double cal_energy() override;
        void cal_force(ModuleBase::matrix& force) override;
        void cal_stress(ModuleBase::matrix& stress) override;
        void after_all_runners() override;

    private:

        // wavefunction coefficients
        // hamilt::Hamilt<T, Device>* p_hamilt = nullptr;

        ModuleDirectMin::DirectMin * dm;


    };  // end class ESolver_DirectMin
} // namespace ModuleESolver

#endif // ESOLVER_DIRECTMIN_H