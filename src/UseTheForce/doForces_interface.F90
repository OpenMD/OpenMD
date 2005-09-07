!! Interfaces for C programs to module....

subroutine initFortranFF(use_RF_c, use_UW_c, use_DW_c, thisStat)
  use doForces, ONLY: init_FF
  integer, intent(in) :: use_RF_c
  integer, intent(in) :: use_UW_c
  integer, intent(in) :: use_DW_c
  integer, intent(out) :: thisStat
  logical :: use_RF, use_UW, use_DW
  
  use_RF = (use_RF_c .ne. 0)
  use_UW = (use_UW_c .ne. 0)
  use_DW = (use_DW_c .ne. 0)

  call init_FF(use_RF, use_UW, use_DW, thisStat)

end subroutine initFortranFF

subroutine doForceloop(q, q_group, A, eFrame, f, t, tau, pot, &
     do_pot_c, do_stress_c, error)

  use definitions, ONLY: dp
  use simulation
  use doForces, ONLY: do_force_loop
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
  real ( kind = dp ) :: pot
  logical ( kind = 2) :: do_pot_c, do_stress_c
  integer :: error

  call do_force_loop(q, q_group, A, eFrame, f, t, tau, pot, &
       do_pot_c, do_stress_c, error)

end subroutine doForceloop

