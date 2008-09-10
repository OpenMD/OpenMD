subroutine setFortranSim(setThisSim, CnGlobal, CnLocal, c_idents, &
     CnExcludes, Cexcludes, CnOneTwo, ConeTwo, CnOneThree, ConeThree, &
     CnOneFour, ConeFour, CmolMembership, Cmfact, CnGroups, &
     CglobalGroupMembership, status)
  use definitions   
  use simulation

  type (simtype) :: setThisSim
  integer, intent(inout) :: CnGlobal, CnLocal
  integer, dimension(CnLocal),intent(inout) :: c_idents

  integer :: CnExcludes
  integer, dimension(2,CnExcludes), intent(inout) :: Cexcludes

  integer :: CnOneTwo
  integer, dimension(2,CnOneTwo), intent(inout) :: ConeTwo

  integer :: CnOneThree
  integer, dimension(2,CnOneThree), intent(inout) :: ConeThree

  integer :: CnOneFour
  integer, dimension(2,CnOneFour), intent(inout) :: ConeFour

  integer, dimension(CnGlobal),intent(inout) :: CmolMembership 
  !!  Result status, success = 0, status = -1
  integer, intent(inout) :: status

  !! mass factors used for molecular cutoffs
  real ( kind = dp ), dimension(CnLocal) :: Cmfact
  integer, intent(in):: CnGroups
  integer, dimension(CnGlobal), intent(inout):: CglobalGroupMembership

  call SimulationSetup(setThisSim, CnGlobal, CnLocal, c_idents, &
       CnExcludes, Cexcludes, CnOneTwo, ConeTwo, CnOneThree, ConeThree, &
       CnOneFour, ConeFour, CmolMembership, Cmfact, CnGroups, &
       CglobalGroupMembership, status)

end subroutine setFortranSim

subroutine setFortranBox(cHmat, cHmatInv, cBoxIsOrthorhombic)
  use simulation, only : setBox
  use definitions
  real(kind=dp), dimension(3,3) :: cHmat, cHmatInv
  integer :: cBoxIsOrthorhombic

  call setBox(cHmat, cHmatInv, cBoxIsOrthorhombic)

end subroutine setFortranBox
