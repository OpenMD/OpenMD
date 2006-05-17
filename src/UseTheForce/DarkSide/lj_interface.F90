subroutine newLJtype(c_ident, sigma, epsilon, soft_pot, status)
  use definitions
  use lj, ONLY : module_newLJtype => newLJtype

  integer,intent(inout) :: c_ident
  real(kind=dp),intent(inout) :: sigma
  real(kind=dp),intent(inout) :: epsilon
  integer,intent(inout) :: soft_pot
  integer,intent(inout) :: status

  call module_newLJtype(c_ident, sigma, epsilon, soft_pot, status)

end subroutine newLJtype

subroutine destroyLJTypes()
  use lj, ONLY: module_destroyLJTypes => destroyLJTypes

  call module_destroyLJTypes()

end subroutine destroyLJTypes
