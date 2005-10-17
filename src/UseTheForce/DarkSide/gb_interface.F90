subroutine newGBtype(c_ident, sigma, l2b_ratio, eps, eps_ratio, mu, nu, &
     status)
  
  use definitions, ONLY : dp
  use gayberne, ONLY : module_newGBtype => newGBtype

  integer, intent(inout) :: c_ident, status
  real( kind = dp ), intent(inout) :: sigma, l2b_ratio, eps, eps_ratio
  real( kind = dp ), intent(inout) :: mu, nu

  call module_newGBtype(c_ident, sigma, l2b_ratio, eps, eps_ratio, &
       mu, nu, status)
  
  return
end subroutine newGBtype

subroutine destroyGBTypes()
  use gayberne, ONLY: module_destroyGBTypes => destroyGBTypes

  call module_destroyGBTypes()

end subroutine destroyGBTypes
