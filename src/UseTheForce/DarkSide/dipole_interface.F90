subroutine newDipoleType(ident, dipole_moment, status)

  use dipole_dipole, ONLY : module_newDipoleType => newDipoleType

  integer, parameter :: DP = selected_real_kind(15)
  integer,intent(inout) :: ident
  real(kind=dp),intent(inout) :: dipole_moment
  integer,intent(inout) :: status
  
  call module_newDipoleType(ident, dipole_moment, status)
  
end subroutine newDipoleType
