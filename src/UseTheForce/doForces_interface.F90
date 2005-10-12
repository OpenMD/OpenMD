!! Interfaces for C programs to module....

subroutine initFortranFF(correctionMethod, thisStat)
  use doForces, ONLY: init_FF
  use definitions, ONLY : dp

  integer, intent(in) :: correctionMethod
  integer, intent(out) :: thisStat
  integer :: correction
  
  correction = correctionMethod
  
  call init_FF(correction, thisStat)

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
  real ( kind = dp ),dimension(POT_ARRAY_SIZE) :: pot
  logical ( kind = 2) :: do_pot_c, do_stress_c
  integer :: error

  call do_force_loop(q, q_group, A, eFrame, f, t, tau, pot, &
       do_pot_c, do_stress_c, error)

end subroutine doForceloop

