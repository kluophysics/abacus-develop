#include "esolver_directmin.h"

namespace ModuleESolver
{

// constructor for DirectMin_LCAO
template <typename TK, typename TR> 
ESolver_DirectMin_LCAO<TK,TR>::ESolver_DirectMin_LCAO()
{
 

}

template <typename TK, typename TR>
ESolver_DirectMin_LCAO<TK,TR>::~ESolver_DirectMin_LCAO()
{

}

template <typename TK, typename TR>
void ESolver_DirectMin_LCAO<TK,TR>::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
  return ;
}

template <typename TK, typename TR>
void ESolver_DirectMin_LCAO<TK,TR>::runner(UnitCell& ucell, const int istep)
{
  return ;
}

template <typename TK, typename TR>
double ESolver_DirectMin_LCAO<TK,TR>::cal_energy( )
{

      return 0;
}

template <typename TK, typename TR>
void ESolver_DirectMin_LCAO<TK,TR>::cal_force(UnitCell& ucell, ModuleBase::matrix& force) 
{

      return;
}

template <typename TK, typename TR>
void ESolver_DirectMin_LCAO<TK,TR>::cal_stress(UnitCell& ucell, ModuleBase::matrix& stress)
{
      return;
}

template <typename TK, typename TR>
void ESolver_DirectMin_LCAO<TK,TR>::after_all_runners(UnitCell& ucell)
{
      return;
}

// template <typename TK, typename TR>
// void ESolver_DirectMin_LCAO<TK,TR>::after_scf(UnitCell& ucell,  const int istep)
// {
//       return;
// }

template <typename TK, typename TR>
void ESolver_DirectMin_LCAO<TK,TR>::others(UnitCell& ucell,  const int istep)
{
      return;
}


    // don't forget to include templates lastly...
    template class ESolver_DirectMin_LCAO<double, double>;
    template class ESolver_DirectMin_LCAO<std::complex<double>, double>;
    template class ESolver_DirectMin_LCAO<std::complex<double>, std::complex<double>>;

}
