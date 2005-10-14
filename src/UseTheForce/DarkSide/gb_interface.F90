subroutine set_gb_pair_params(sigma, l2b_ratio, eps, eps_ratio, mu, nu)

  use definitions, ONLY : dp
  use gb_pair, ONLY : module_set_gb_pair_params => set_gb_pair_params

  real( kind = dp ), intent(inout) :: sigma, l2b_ratio, eps, eps_ratio
  real( kind = dp ), intent(inout) :: mu, nu

  call module_set_gb_pair_params(sigma, l2b_ratio, eps, eps_ratio, mu, nu)

  return
end subroutine set_gb_pair_params

subroutine completeGayBerneFF(status)
  use gb_pair, only: complete_GayBerne_FF
  
  intent(out) :: status
  integer :: myStatus

  myStatus = 0
  
  call complete_GayBerne_FF(myStatus)

  status = myStatus

  return
end subroutine completeGayBerneFF

subroutine destroyGayBerneTypes()

  use gb_pair, only: module_destroyGayBerneTypes => destroyGayBerneTypes
  call module_destroyGayBerneTypes()
 
end subroutine destroyGayBerneTypes
