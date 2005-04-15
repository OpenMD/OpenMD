subroutine notifyFortranCutoffs(this_rcut, this_rsw, this_rlist )
  use notifyCutoffs
  use definitions, ONLY : dp

  real(kind=dp), intent(in) :: this_rcut, this_rsw, this_rlist

  call cutoffNotify(this_rcut, this_rsw, this_rlist )

end subroutine notifyFortranCutoffs
