subroutine newLJtype(c_ident, sigma, epsilon, status)
  use lj, ONLY : module_newLJtype => newLJtype
  integer, parameter :: DP = selected_real_kind(15)
  integer,intent(inout) :: c_ident
  real(kind=dp),intent(inout) :: sigma
  real(kind=dp),intent(inout) :: epsilon
  integer,intent(inout) :: status
  
  call module_newLJtype(c_ident, sigma, epsilon, status)
  
end subroutine newLJtype

subroutine useGeometricMixing()
  use lj, ONLY: module_useGeometricMixing => useGeometricMixing
  
  call module_useGeometricMixing()
  return
end subroutine useGeometricMixing
