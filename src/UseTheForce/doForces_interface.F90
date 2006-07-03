!! Interfaces for C programs to module....

subroutine initFortranFF(thisStat)
  use doForces, ONLY: init_FF
  use definitions, ONLY : dp

  integer, intent(out) :: thisStat
  
  call init_FF(thisStat)

end subroutine initFortranFF

subroutine doForceloop(q, q_group, A, eFrame, f, t, tau, pot, &
     do_pot_c, do_stress_c, error)

  use definitions, ONLY: dp
  use simulation
  use doForces, ONLY: do_force_loop

#define __FORTRAN90
#include "UseTheForce/DarkSide/fInteractionMap.h"

  !! Position array provided by C, dimensioned by getNlocal
  real ( kind = dp ), dimension(3, nLocal) :: q
  !! molecular center-of-mass position array
  real ( kind = dp ), dimension(3, nGroups) :: q_group
  !! Rotation Matrix for each long range particle in simulation.
  real( kind = dp), dimension(9, nLocal) :: A    
  !! Unit vectors for dipoles (lab frame)
  real( kind = dp ), dimension(9,nLocal) :: eFrame
  !! Force array provided by C, dimensioned by getNlocal
  real ( kind = dp ), dimension(3,nLocal) :: f
  !! Torsion array provided by C, dimensioned by getNlocal
  real( kind = dp ), dimension(3,nLocal) :: t    

  !! Stress Tensor
  real( kind = dp), dimension(9) :: tau   
  real ( kind = dp ),dimension(LR_POT_TYPES) :: pot
  logical ( kind = 2) :: do_pot_c, do_stress_c
  integer :: error

  call do_force_loop(q, q_group, A, eFrame, f, t, tau, pot, &
       do_pot_c, do_stress_c, error)

end subroutine doForceloop

subroutine getAccumulatedBoxDipole( box_dipole )
  
  use definitions, ONLY: dp
  use doForces, ONLY: getBoxDipole

  !! simulation box dipole moment 
  real ( kind = dp ), dimension(3) :: box_dipole
  
  call getBoxDipole( box_dipole )
  
end subroutine getAccumulatedBoxDipole

subroutine setAccumulateBoxDipole()

  use doForces, ONLY: setBoxDipole

  call setBoxDipole()
  
end subroutine setAccumulateBoxDipole

subroutine setFortranElectrostaticMethod(electrostaticMethod)
  use doForces, ONLY : setElectrostaticMethod

  integer, intent(in) :: electrostaticMethod

  call setElectrostaticMethod(electrostaticMethod)

end subroutine setFortranElectrostaticMethod

subroutine notifyFortranCutoffPolicy(cutPolicy)
  use doForces, ONLY : setCutoffPolicy

  integer, intent(in) :: cutPolicy

  call setCutoffPolicy( cutPolicy )

end subroutine notifyFortranCutoffPolicy

subroutine notifyFortranSkinThickness(this_skin)
  use doForces, ONLY : setSkinThickness
  use definitions, ONLY : dp

  real(kind=dp), intent(in) :: this_skin

  call setSkinThickness( this_skin )

end subroutine notifyFortranSkinThickness

subroutine notifyFortranCutoffs(this_rcut, this_rsw)
  use doForces, ONLY : setCutoffs
  use definitions, ONLY : dp

  real(kind=dp), intent(in) :: this_rcut, this_rsw

  call setCutoffs(this_rcut, this_rsw)

end subroutine notifyFortranCutoffs

subroutine notifyFortranYouAreOnYourOwn()
  use doForces, ONLY : cWasLame

  call cWasLame()
end subroutine notifyFortranYouAreOnYourOwn
