subroutine set_gb_pair_params(sigma, l2b_ratio, eps, eps_ratio, mu, nu)

  use definitions, ONLY : dp
  use gb_pair, ONLY : module_set_gb_pair_params => set_gb_pair_params

  real( kind = dp ), intent(inout) :: sigma, l2b_ratio, eps, eps_ratio
  real( kind = dp ), intent(inout) :: mu, nu

  call module_set_gb_pair_params(sigma, l2b_ratio, eps, eps_ratio, mu, nu)

end subroutine set_gb_pair_params
