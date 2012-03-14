subroutine newGBtype(c_ident, d, l, epsX, epsS, epsE, dw, status)
  
  use definitions, ONLY : dp
  use gayberne, ONLY : module_newGBtype => newGBtype

  integer, intent(inout) :: c_ident, status
  real( kind = dp ), intent(inout) :: d, l, epsX, epsS, epsE, dw

  call module_newGBtype(c_ident, d, l, epsX, epsS, epsE, dw, status)
  
  return
end subroutine newGBtype

subroutine completeGBFF(status)

  use gayberne, only: complete_GB_FF

  integer, intent(out)  :: status
  integer :: myStatus

  myStatus = 0

  call complete_GB_FF(myStatus)

  status = myStatus

  return
end subroutine completeGBFF


subroutine destroyGBTypes()
  use gayberne, ONLY: module_destroyGBTypes => destroyGBTypes

  call module_destroyGBTypes()

end subroutine destroyGBTypes
