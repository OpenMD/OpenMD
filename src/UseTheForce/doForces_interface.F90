!! Interfaces for C programs to module....

subroutine initFortranFF(thisStat)
  use doForces, ONLY: init_FF
  use definitions, ONLY : dp

  integer, intent(out) :: thisStat
  
  call init_FF(thisStat)

end subroutine initFortranFF

subroutine doForceloop(q, q_group, A, eFrame, f, t, tau, pot, particle_pot, &
     error)
  
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
  real( kind = dp ), dimension(nLocal) :: particle_pot
  integer :: error

  call do_force_loop(q, q_group, A, eFrame, f, t, tau, pot, particle_pot, &
       error)
  
end subroutine doForceloop

subroutine notifyFortranSkinThickness(this_skin)
  use doForces, ONLY : setSkinThickness
  use definitions, ONLY : dp

  real(kind=dp), intent(in) :: this_skin

  call setSkinThickness( this_skin )

end subroutine notifyFortranSkinThickness

subroutine notifyFortranCutoffs(this_rcut, this_rsw, this_sp, this_sf)
  use doForces, ONLY : setCutoffs
  use definitions, ONLY : dp

  real(kind=dp), intent(in) :: this_rcut, this_rsw
  integer, intent(in) :: this_sp, this_sf

  call setCutoffs(this_rcut, this_rsw, this_sp, this_sf)

end subroutine notifyFortranCutoffs
