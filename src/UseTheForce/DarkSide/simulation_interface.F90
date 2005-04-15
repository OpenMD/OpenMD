subroutine setFortranSim(setThisSim, CnGlobal, CnLocal, c_idents, &
     CnLocalExcludes, CexcludesLocal, CnGlobalExcludes, CexcludesGlobal, &
     CmolMembership, Cmfact, CnGroups, CglobalGroupMembership, &
     status)
  use definitions, ONLY : dp    
  use simulation

  type (simtype) :: setThisSim
  integer, intent(inout) :: CnGlobal, CnLocal
  integer, dimension(CnLocal),intent(inout) :: c_idents

  integer :: CnLocalExcludes
  integer, dimension(2,CnLocalExcludes), intent(inout) :: CexcludesLocal
  integer :: CnGlobalExcludes
  integer, dimension(CnGlobalExcludes), intent(inout) :: CexcludesGlobal
  integer, dimension(CnGlobal),intent(inout) :: CmolMembership 
  !!  Result status, success = 0, status = -1
  integer, intent(inout) :: status

  !! mass factors used for molecular cutoffs
  real ( kind = dp ), dimension(CnLocal) :: Cmfact
  integer, intent(in):: CnGroups
  integer, dimension(CnGlobal), intent(inout):: CglobalGroupMembership

  call SimulationSetup(setThisSim, CnGlobal, CnLocal, c_idents, &
       CnLocalExcludes, CexcludesLocal, CnGlobalExcludes, CexcludesGlobal, &
       CmolMembership, Cmfact, CnGroups, CglobalGroupMembership, &
       status)
end subroutine setFortranSim

subroutine setFortranBox(cHmat, cHmatInv, cBoxIsOrthorhombic)
  use simulation, only : setBox
  use definitions, ONLY : dp
  real(kind=dp), dimension(3,3) :: cHmat, cHmatInv
  integer :: cBoxIsOrthorhombic

  call setBox(cHmat, cHmatInv, cBoxIsOrthorhombic)

end subroutine setFortranBox
