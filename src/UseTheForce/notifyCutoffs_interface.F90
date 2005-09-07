subroutine notifyFortranCutoffs(this_rcut, this_rsw, this_rlist, cutPolicy)
  use notifyCutoffs
  use definitions, ONLY : dp

  real(kind=dp), intent(in) :: this_rcut, this_rsw, this_rlist
  integer, intent(in) :: cutPolicy

  call cutoffNotify(this_rcut, this_rsw, this_rlist, cutPolicy )

end subroutine notifyFortranCutoffs
