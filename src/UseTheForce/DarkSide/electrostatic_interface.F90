subroutine setElectrostaticSummationMethod(the_ESM)
  use electrostatic_module, ONLY : module_setESM => setElectrostaticSummationMethod
  integer,intent(inout) :: the_ESM
  call module_setESM(the_ESM)
end subroutine setElectrostaticSummationMethod

subroutine setScreeningMethod(the_SM)
  use electrostatic_module, ONLY : module_setSM => setScreeningMethod
  integer,intent(inout) :: the_SM
  call module_setSM(the_SM)
end subroutine setScreeningMethod

subroutine setElectrostaticCutoffRadius(the_rcut, the_rsw)

  use definitions, ONLY : dp
  use electrostatic_module, ONLY : module_setECR => setElectrostaticCutoffRadius

  real(kind=dp), intent(inout) :: the_rcut
  real(kind=dp), intent(inout) :: the_rsw
  call module_setECR(the_rcut, the_rsw)

end subroutine setElectrostaticCutoffRadius

subroutine setDampingAlpha(the_alpha)

  use definitions, ONLY : dp
  use electrostatic_module, ONLY : module_setDA => setDampingAlpha

  real(kind=dp),intent(inout) :: the_alpha
  call module_setDA(the_alpha)

end subroutine setDampingAlpha
 
subroutine setReactionFieldDielectric(the_dielectric)

  use definitions, ONLY : dp
  use electrostatic_module, ONLY : module_setRFD => setReactionFieldDielectric

  real(kind=dp),intent(inout) :: the_dielectric
  call module_setRFD(the_dielectric)

end subroutine setReactionFieldDielectric

subroutine newElectrostaticType(atp, status)

  use electrostatic_module, ONLY : module_newElectrostaticType => newElectrostaticType

#define __FORTRAN90
#include "types/AtomTypeProperties.h"

  type(AtomTypeProperties), intent(in) :: atp
  integer, intent(inout) :: status

  integer :: ident
  logical :: is_Electrostatic, is_Charge, is_Dipole
  logical :: is_SplitDipole, is_Quadrupole, is_Tap

  ident = atp%ident
  is_Electrostatic = ((atp%is_Charge .ne. 0) .or. &
       (atp%is_Dipole .ne. 0)) .or. &
       (atp%is_Quadrupole .ne. 0)
  is_Charge = (atp%is_Charge .ne. 0)
  is_Dipole = (atp%is_Dipole .ne. 0)
  is_SplitDipole = (atp%is_SplitDipole .ne. 0)
  is_Quadrupole = (atp%is_Quadrupole .ne. 0)
  is_Tap = (atp%is_StickyPower .ne. 0)

  call module_newElectrostaticType(ident, is_Charge, is_Dipole, &
       is_SplitDipole, is_Quadrupole, is_Tap, status)

end subroutine newElectrostaticType

subroutine setCharge(ident, charge, status)

  use definitions, ONLY : dp
  use electrostatic_module, ONLY : module_setCharge => setCharge

  integer,intent(inout) :: ident
  real(kind=dp),intent(inout) :: charge
  integer,intent(inout) :: status

  call module_setCharge(ident, charge, status)

end subroutine setCharge

subroutine setDipoleMoment(ident, dipole_moment, status)

  use definitions, ONLY : dp
  use electrostatic_module, ONLY : module_setDipoleMoment => setDipoleMoment

  integer,intent(inout) :: ident
  real(kind=dp),intent(inout) :: dipole_moment
  integer,intent(inout) :: status

  call module_setDipoleMoment(ident, dipole_moment, status)

end subroutine setDipoleMoment

subroutine setSplitDipoleDistance(ident, split_dipole_distance, status)

  use definitions, ONLY : dp
  use electrostatic_module, ONLY : module_setSplitDipoleDistance => setSplitDipoleDistance

  integer,intent(inout) :: ident
  real(kind=dp),intent(inout) :: split_dipole_distance
  integer,intent(inout) :: status

  call module_setSplitDipoleDistance(ident, split_dipole_distance, status)

end subroutine setSplitDipoleDistance

subroutine setQuadrupoleMoments(ident, quadrupole_moments, status)

  use definitions, ONLY : dp
  use electrostatic_module, ONLY : module_setQuadrupoleMoments => setQuadrupoleMoments

  integer,intent(inout) :: ident
  real(kind=dp),intent(inout),dimension(3) :: quadrupole_moments
  integer,intent(inout) :: status

  call module_setQuadrupoleMoments(ident, quadrupole_moments, status)

end subroutine setQuadrupoleMoments

subroutine destroyElectrostaticTypes()
  use electrostatic_module, ONLY: m_destroyElectrostaticTypes =>destroyElectrostaticTypes

  call m_destroyElectrostaticTypes()

end subroutine destroyElectrostaticTypes
