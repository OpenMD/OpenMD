!!
!! Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
!!
!! The University of Notre Dame grants you ("Licensee") a
!! non-exclusive, royalty free, license to use, modify and
!! redistribute this software in source and binary code form, provided
!! that the following conditions are met:
!!
!! 1. Redistributions of source code must retain the above copyright
!!    notice, this list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above copyright
!!    notice, this list of conditions and the following disclaimer in the
!!    documentation and/or other materials provided with the
!!    distribution.
!!
!! This software is provided "AS IS," without a warranty of any
!! kind. All express or implied conditions, representations and
!! warranties, including any implied warranty of merchantability,
!! fitness for a particular purpose or non-infringement, are hereby
!! excluded.  The University of Notre Dame and its licensors shall not
!! be liable for any damages suffered by licensee as a result of
!! using, modifying or distributing the software or its
!! derivatives. In no event will the University of Notre Dame or its
!! licensors be liable for any lost revenue, profit or data, or for
!! direct, indirect, special, consequential, incidental or punitive
!! damages, however caused and regardless of the theory of liability,
!! arising out of the use of or inability to use software, even if the
!! University of Notre Dame has been advised of the possibility of
!! such damages.
!!
!! SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
!! research, please cite the appropriate papers when you publish your
!! work.  Good starting points are:
!!                                                                      
!! [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
!! [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
!! [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
!! [4]  Vardeman & Gezelter, in progress (2009).
!!

!! Fortran interface to C entry plug.

module simulation
  use definitions
  use status
  use linearAlgebra
  use neighborLists
  use force_globals
  use vector_class
  use atype_module
  use switcheroo
#ifdef IS_MPI
  use mpiSimulation
#endif

  implicit none
  PRIVATE

#define __FORTRAN90
#include "brains/fSimulation.h"
#include "UseTheForce/fSwitchingFunction.h"

  type (simtype), public, save :: thisSim

  logical, save :: simulation_setup_complete = .false.

  integer, public, save :: nLocal, nGlobal
  integer, public, save :: nGroups, nGroupGlobal
  integer, public, save :: nExcludes = 0
  integer, public, save :: nOneTwo = 0
  integer, public, save :: nOneThree = 0
  integer, public, save :: nOneFour = 0

  integer, allocatable, dimension(:,:), public :: excludes
  integer, allocatable, dimension(:),   public :: molMembershipList
  integer, allocatable, dimension(:),   public :: groupListRow
  integer, allocatable, dimension(:),   public :: groupStartRow
  integer, allocatable, dimension(:),   public :: groupListCol
  integer, allocatable, dimension(:),   public :: groupStartCol
  integer, allocatable, dimension(:),   public :: groupListLocal
  integer, allocatable, dimension(:),   public :: groupStartLocal
  integer, allocatable, dimension(:),   public :: nSkipsForLocalAtom
  integer, allocatable, dimension(:,:), public :: skipsForLocalAtom
  integer, allocatable, dimension(:),   public :: nSkipsForRowAtom
  integer, allocatable, dimension(:,:), public :: skipsForRowAtom
  integer, allocatable, dimension(:),   public :: nTopoPairsForAtom
  integer, allocatable, dimension(:,:),   public :: toposForAtom
  integer, allocatable, dimension(:,:),   public :: topoDistance
  real(kind=dp), allocatable, dimension(:), public :: mfactRow
  real(kind=dp), allocatable, dimension(:), public :: mfactCol
  real(kind=dp), allocatable, dimension(:), public :: mfactLocal

  real(kind=dp), public, dimension(3,3), save :: Hmat, HmatInv
  real(kind=dp), save :: DangerRcut
  logical, public, save :: boxIsOrthorhombic

  public :: SimulationSetup
  public :: getNlocal
  public :: setBox
  public :: checkBox
  public :: SimUsesPBC
  public :: SimUsesAtomicVirial
  public :: SimUsesDirectionalAtoms
  public :: SimUsesMetallicAtoms
  public :: SimRequiresSkipCorrection
  public :: SimRequiresSelfCorrection
  public :: setHmatDangerousRcutValue

contains

  subroutine SimulationSetup(setThisSim, CnGlobal, CnLocal, c_idents, &
       CnExcludes, Cexcludes, CnOneTwo, ConeTwo, CnOneThree, ConeThree, &
       CnOneFour, ConeFour, CmolMembership, Cmfact, CnGroups, &
       CglobalGroupMembership, status)    

    type (simtype) :: setThisSim
    integer, intent(inout) :: CnGlobal, CnLocal
    integer, dimension(CnLocal), intent(inout) :: c_idents

    integer :: CnExcludes
    integer, dimension(2,CnExcludes), intent(in) :: Cexcludes
    integer :: CnOneTwo
    integer, dimension(2,CnOneTwo), intent(in) :: ConeTwo
    integer :: CnOneThree
    integer, dimension(2,CnOneThree), intent(in) :: ConeThree
    integer :: CnOneFour
    integer, dimension(2,CnOneFour), intent(in) :: ConeFour

    integer, dimension(CnGlobal),intent(in) :: CmolMembership 
    !!  Result status, success = 0, status = -1
    integer, intent(out) :: status
    integer :: i, j, me, thisStat, alloc_stat, myNode, id1, id2
    integer :: ia, jend

    !! mass factors used for molecular cutoffs
    real ( kind = dp ), dimension(CnLocal) :: Cmfact
    integer, intent(in):: CnGroups
    integer, dimension(CnGlobal), intent(in):: CglobalGroupMembership
    integer :: maxSkipsForLocalAtom, maxToposForAtom, glPointer
    integer :: maxSkipsForRowAtom

#ifdef IS_MPI
    integer, allocatable, dimension(:) :: c_idents_Row
    integer, allocatable, dimension(:) :: c_idents_Col
    integer :: nAtomsInRow, nGroupsInRow, aid
    integer :: nAtomsInCol, nGroupsInCol, gid
#endif  

    simulation_setup_complete = .false.
    status = 0

    ! copy C struct into fortran type

    nLocal = CnLocal
    nGlobal = CnGlobal
    nGroups = CnGroups

    thisSim = setThisSim

    nExcludes = CnExcludes

    call InitializeForceGlobals(nLocal, thisStat)
    if (thisStat /= 0) then
       write(default_error,*) "SimSetup: InitializeForceGlobals error"
       status = -1
       return
    endif

    call InitializeSimGlobals(thisStat)
    if (thisStat /= 0) then
       write(default_error,*) "SimSetup: InitializeSimGlobals error"
       status = -1
       return
    endif

#ifdef IS_MPI
    ! We can only set up forces if mpiSimulation has been setup.
    if (.not. isMPISimSet()) then
       write(default_error,*) "MPI is not set"
       status = -1
       return
    endif
    nAtomsInRow = getNatomsInRow(plan_atom_row)
    nAtomsInCol = getNatomsInCol(plan_atom_col)
    nGroupsInRow = getNgroupsInRow(plan_group_row)
    nGroupsInCol = getNgroupsInCol(plan_group_col)
    mynode = getMyNode()

    call gather(c_idents, c_idents_Row, plan_atom_row)
    call gather(c_idents, c_idents_Col, plan_atom_col)

    do i = 1, nAtomsInRow
       me = getFirstMatchingElement(atypes, "c_ident", c_idents_Row(i))
       atid_Row(i) = me
    enddo

    do i = 1, nAtomsInCol
       me = getFirstMatchingElement(atypes, "c_ident", c_idents_Col(i))
       atid_Col(i) = me
    enddo

#endif

#ifdef IS_MPI
    allocate(groupStartRow(nGroupsInRow+1),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1
       return
    endif
    allocate(groupStartCol(nGroupsInCol+1),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1
       return
    endif
    allocate(groupListRow(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1
       return
    endif
    allocate(groupListCol(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1
       return
    endif
    allocate(mfactRow(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1
       return
    endif
    allocate(mfactCol(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1
       return
    endif
    allocate(mfactLocal(nLocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1
       return
    endif
    
    glPointer = 1
    
    do i = 1, nGroupsInRow 
       
       gid = GroupRowToGlobal(i)
       groupStartRow(i) = glPointer       
       
       do j = 1, nAtomsInRow
          aid = AtomRowToGlobal(j)
          if (CglobalGroupMembership(aid) .eq. gid) then
             groupListRow(glPointer) = j
             glPointer = glPointer + 1
          endif
       enddo
    enddo
    groupStartRow(nGroupsInRow+1) = nAtomsInRow + 1
    
    glPointer = 1
    
    do i = 1, nGroupsInCol
       
       gid = GroupColToGlobal(i)
       groupStartCol(i) = glPointer       
       
       do j = 1, nAtomsInCol
          aid = AtomColToGlobal(j)
          if (CglobalGroupMembership(aid) .eq. gid) then
             groupListCol(glPointer) = j
             glPointer = glPointer + 1
          endif
       enddo
    enddo
    groupStartCol(nGroupsInCol+1) = nAtomsInCol + 1

    mfactLocal = Cmfact        

    call gather(mfactLocal,      mfactRow,      plan_atom_row)
    call gather(mfactLocal,      mfactCol,      plan_atom_col)

    if (allocated(mfactLocal)) then
       deallocate(mfactLocal)
    end if
#else
    allocate(groupStartRow(nGroups+1),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1
       return
    endif
    allocate(groupStartCol(nGroups+1),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1
       return
    endif
    allocate(groupListRow(nLocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1
       return
    endif
    allocate(groupListCol(nLocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1
       return
    endif
    allocate(mfactRow(nLocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1
       return
    endif
    allocate(mfactCol(nLocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1
       return
    endif
    allocate(mfactLocal(nLocal),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1
       return
    endif

    glPointer = 1
    do i = 1, nGroups
       groupStartRow(i) = glPointer       
       groupStartCol(i) = glPointer
       do j = 1, nLocal
          if (CglobalGroupMembership(j) .eq. i) then
             groupListRow(glPointer) = j
             groupListCol(glPointer) = j
             glPointer = glPointer + 1
          endif
       enddo
    enddo
    groupStartRow(nGroups+1) = nLocal + 1
    groupStartCol(nGroups+1) = nLocal + 1

    do i = 1, nLocal
       mfactRow(i) = Cmfact(i)
       mfactCol(i) = Cmfact(i)
    end do

#endif

    ! We build the local atid's for both mpi and nonmpi
    do i = 1, nLocal
       me = getFirstMatchingElement(atypes, "c_ident", c_idents(i))
       atid(i) = me
       c_idents_local(i) = c_idents(i)
    enddo

    do i = 1, nExcludes
       excludes(1,i) = Cexcludes(1,i)
       excludes(2,i) = Cexcludes(2,i)
    enddo

#ifdef IS_MPI
    allocate(nSkipsForRowAtom(nAtomsInRow), stat=alloc_stat)
#endif

    allocate(nSkipsForLocalAtom(nLocal), stat=alloc_stat)
    
    if (alloc_stat /= 0 ) then
       thisStat = -1
       write(*,*) 'Could not allocate nSkipsForAtom array'
       return
    endif

#ifdef IS_MPI
    maxSkipsForRowAtom = 0
    do j = 1, nAtomsInRow
       nSkipsForRowAtom(j) = 0
       id1 = AtomRowToGlobal(j)
       do i = 1, nExcludes
          if (excludes(1,i) .eq. id1 ) then
             nSkipsForRowAtom(j) = nSkipsForRowAtom(j) + 1             
             if (nSkipsForRowAtom(j) .gt. maxSkipsForRowAtom) then
                maxSkipsForRowAtom = nSkipsForRowAtom(j)
             endif
          endif
          if (excludes(2,i) .eq. id1 ) then
             nSkipsForRowAtom(j) = nSkipsForRowAtom(j) + 1             
             if (nSkipsForRowAtom(j) .gt. maxSkipsForRowAtom) then
                maxSkipsForRowAtom = nSkipsForRowAtom(j)
             endif
          endif
       end do
    enddo
#endif
    maxSkipsForLocalAtom = 0
    do j = 1, nLocal
       nSkipsForLocalAtom(j) = 0
#ifdef IS_MPI
       id1 = AtomLocalToGlobal(j)
#else
       id1 = j
#endif
       do i = 1, nExcludes
          if (excludes(1,i) .eq. id1 ) then
             nSkipsForLocalAtom(j) = nSkipsForLocalAtom(j) + 1             
             if (nSkipsForLocalAtom(j) .gt. maxSkipsForLocalAtom) then
                maxSkipsForLocalAtom = nSkipsForLocalAtom(j)
             endif
          endif
          if (excludes(2,i) .eq. id1 ) then
             nSkipsForLocalAtom(j) = nSkipsForLocalAtom(j) + 1             
             if (nSkipsForLocalAtom(j) .gt. maxSkipsForLocalAtom) then
                maxSkipsForLocalAtom = nSkipsForLocalAtom(j)
             endif
          endif
       end do
    enddo
    
#ifdef IS_MPI
    allocate(skipsForRowAtom(nAtomsInRow, maxSkipsForRowAtom), stat=alloc_stat)
#endif
    allocate(skipsForLocalAtom(nLocal, maxSkipsForLocalAtom), stat=alloc_stat)

    if (alloc_stat /= 0 ) then
       write(*,*) 'Could not allocate skipsForAtom arrays'
       return
    endif
    
#ifdef IS_MPI
    do j = 1, nAtomsInRow
       nSkipsForRowAtom(j) = 0 
       id1 = AtomRowToGlobal(j)
       do i = 1, nExcludes
          if (excludes(1,i) .eq. id1 ) then
             nSkipsForRowAtom(j) = nSkipsForRowAtom(j) + 1
             ! exclude lists have global ID's 
             id2 = excludes(2,i)
             skipsForRowAtom(j, nSkipsForRowAtom(j)) = id2
          endif
          if (excludes(2, i) .eq. id1 ) then
             nSkipsForRowAtom(j) = nSkipsForRowAtom(j) + 1
             ! exclude lists have global ID's 
             id2 = excludes(1,i)
             skipsForRowAtom(j, nSkipsForRowAtom(j)) = id2
          endif
       end do
    enddo
#endif
    do j = 1, nLocal
       nSkipsForLocalAtom(j) = 0 
#ifdef IS_MPI
       id1 = AtomLocalToGlobal(j)
#else
       id1 = j
#endif
       do i = 1, nExcludes
          if (excludes(1,i) .eq. id1 ) then
             nSkipsForLocalAtom(j) = nSkipsForLocalAtom(j) + 1
             ! exclude lists have global ID's
#ifdef IS_MPI
             id2 = AtomGlobalToLocal(excludes(2,i))
#else
             id2 = excludes(2,i)
#endif
             skipsForLocalAtom(j, nSkipsForLocalAtom(j)) = id2
          endif
          if (excludes(2, i) .eq. id1 ) then
             nSkipsForLocalAtom(j) = nSkipsForLocalAtom(j) + 1
             ! exclude lists have global ID's 
#ifdef IS_MPI
             id2 = AtomGlobalToLocal(excludes(1,i))
#else
             id2 = excludes(1,i)
#endif
             skipsForLocalAtom(j, nSkipsForLocalAtom(j)) = id2
          endif
       end do
    enddo
           
    do i = 1, nGlobal
       molMemberShipList(i) = CmolMembership(i)
    enddo

#ifdef IS_MPI
    allocate(nTopoPairsForAtom(nAtomsInRow), stat=alloc_stat)
#else
    allocate(nTopoPairsForAtom(nLocal), stat=alloc_stat)
#endif
    if (alloc_stat /= 0 ) then
       thisStat = -1
       write(*,*) 'Could not allocate nTopoPairsForAtom array'
       return
    endif

#ifdef IS_MPI
    jend = nAtomsInRow
#else
    jend = nLocal
#endif
    
    do j = 1, jend
       nTopoPairsForAtom(j) = 0
#ifdef IS_MPI
       id1 = AtomRowToGlobal(j)
#else 
       id1 = j
#endif
       do i = 1, CnOneTwo
          if (ConeTwo(1,i) .eq. id1 ) then
             nTopoPairsForAtom(j) = nTopoPairsForAtom(j) + 1
          endif
          if (ConeTwo(2,i) .eq. id1 ) then
             nTopoPairsForAtom(j) = nTopoPairsForAtom(j) + 1
          endif
       end do

       do i = 1, CnOneThree
          if (ConeThree(1,i) .eq. id1 ) then
             nTopoPairsForAtom(j) = nTopoPairsForAtom(j) + 1
          endif
          if (ConeThree(2,i) .eq. id1 ) then
             nTopoPairsForAtom(j) = nTopoPairsForAtom(j) + 1
          endif
       end do

       do i = 1, CnOneFour
          if (ConeFour(1,i) .eq. id1 ) then
             nTopoPairsForAtom(j) = nTopoPairsForAtom(j) + 1
          endif
          if (ConeFour(2,i) .eq. id1 ) then
             nTopoPairsForAtom(j) = nTopoPairsForAtom(j) + 1
          endif
       end do
    enddo
    
    maxToposForAtom = maxval(nTopoPairsForAtom)
#ifdef IS_MPI
    allocate(toposForAtom(nAtomsInRow, maxToposForAtom), stat=alloc_stat)
    allocate(topoDistance(nAtomsInRow, maxToposForAtom), stat=alloc_stat)
#else
    allocate(toposForAtom(nLocal, maxToposForAtom), stat=alloc_stat)
    allocate(topoDistance(nLocal, maxToposForAtom), stat=alloc_stat)
#endif
    if (alloc_stat /= 0 ) then
       write(*,*) 'Could not allocate topoDistance array'
       return
    endif
    
#ifdef IS_MPI
    jend = nAtomsInRow
#else
    jend = nLocal
#endif
    do j = 1, jend
       nTopoPairsForAtom(j) = 0
#ifdef IS_MPI
       id1 = AtomRowToGlobal(j)
#else 
       id1 = j
#endif
       do i = 1, CnOneTwo
          if (ConeTwo(1,i) .eq. id1 ) then
             nTopoPairsForAtom(j) = nTopoPairsForAtom(j) + 1
             id2 = ConeTwo(2,i)
             toposForAtom(j, nTopoPairsForAtom(j)) = id2
             topoDistance(j, nTopoPairsForAtom(j)) = 1
          endif
          if (ConeTwo(2, i) .eq. id1 ) then
             nTopoPairsForAtom(j) = nTopoPairsForAtom(j) + 1
             id2 = ConeTwo(1,i)
             toposForAtom(j, nTopoPairsForAtom(j)) = id2
             topoDistance(j, nTopoPairsForAtom(j)) = 1
          endif
       end do

       do i = 1, CnOneThree
          if (ConeThree(1,i) .eq. id1 ) then
             nTopoPairsForAtom(j) = nTopoPairsForAtom(j) + 1
             id2 = ConeThree(2,i)
             toposForAtom(j, nTopoPairsForAtom(j)) = id2
             topoDistance(j, nTopoPairsForAtom(j)) = 2
          endif
          if (ConeThree(2, i) .eq. id1 ) then
             nTopoPairsForAtom(j) = nTopoPairsForAtom(j) + 1
             id2 = ConeThree(1,i)
             toposForAtom(j, nTopoPairsForAtom(j)) = id2
             topoDistance(j, nTopoPairsForAtom(j)) = 2
          endif
       end do

       do i = 1, CnOneFour
          if (ConeFour(1,i) .eq. id1 ) then
             nTopoPairsForAtom(j) = nTopoPairsForAtom(j) + 1
             id2 = ConeFour(2,i)
             toposForAtom(j, nTopoPairsForAtom(j)) = id2
             topoDistance(j, nTopoPairsForAtom(j)) = 3
          endif
          if (ConeFour(2, i) .eq. id1 ) then
             nTopoPairsForAtom(j) = nTopoPairsForAtom(j) + 1
             id2 = ConeFour(1,i)
             toposForAtom(j, nTopoPairsForAtom(j)) = id2
             topoDistance(j, nTopoPairsForAtom(j)) = 3
          endif
       end do
    enddo
    
    if (status == 0) simulation_setup_complete = .true.
    
  end subroutine SimulationSetup

  subroutine setBox(cHmat, cHmatInv, cBoxIsOrthorhombic)
    real(kind=dp), dimension(3,3) :: cHmat, cHmatInv
    integer :: cBoxIsOrthorhombic
    integer :: smallest, status
    
    Hmat = cHmat
    HmatInv = cHmatInv
    if (cBoxIsOrthorhombic .eq. 0 ) then
       boxIsOrthorhombic = .false.
    else 
       boxIsOrthorhombic = .true.
    endif
    
    call checkBox()
    return
  end subroutine setBox
  
        subroutine checkBox()

          integer :: i
          real(kind=dp), dimension(3) :: hx, hy, hz, ax, ay, az, piped
          character(len = statusMsgSize) :: errMsg

          hx = Hmat(1,:)
          hy = Hmat(2,:)
          hz = Hmat(3,:)

          ax = cross_product(hy, hz)
          ay = cross_product(hx, hz)
          az = cross_product(hx, hy)

          ax = ax / length(ax)
          ay = ay / length(ay)
          az = az / length(az)

          piped(1) = abs(dot_product(ax, hx))
          piped(2) = abs(dot_product(ay, hy))
          piped(3) = abs(dot_product(az, hz))
          
          do i = 1, 3
             if ((0.5_dp * piped(i)).lt.DangerRcut) then                
                write(errMsg, '(a94,f9.4,a1)') 'One of the dimensions of the Periodic ' &
                     // 'Box is smaller than ' // newline // tab // &
                     'the largest cutoff radius' // &
                     ' (rCut = ', DangerRcut, ')'
                call handleError("checkBox", errMsg)
				
             end if
          enddo        
          return     
        end subroutine checkBox
        
        function SimUsesPBC() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_PBC
        end function SimUsesPBC
        
        function SimUsesAtomicVirial() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_AtomicVirial
        end function SimUsesAtomicVirial
        
        function SimUsesDirectionalAtoms() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_DirectionalAtoms
        end function SimUsesDirectionalAtoms
        
        function SimUsesMetallicAtoms() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_MetallicAtoms
        end function SimUsesMetallicAtoms
        
        function SimRequiresSkipCorrection() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_requires_SkipCorrection
        end function SimRequiresSkipCorrection
        
        function SimRequiresSelfCorrection() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_requires_SelfCorrection
        end function SimRequiresSelfCorrection
        
        subroutine InitializeSimGlobals(thisStat)
          integer, intent(out) :: thisStat
          integer :: alloc_stat
          
          thisStat = 0
          
          call FreeSimGlobals()    
          
          allocate(excludes(2,nExcludes), stat=alloc_stat)
          if (alloc_stat /= 0 ) then
             thisStat = -1
             return
          endif

          allocate(molMembershipList(nGlobal), stat=alloc_stat)
          if (alloc_stat /= 0 ) then
             thisStat = -1
             return
          endif

        end subroutine InitializeSimGlobals

        subroutine FreeSimGlobals()

          !We free in the opposite order in which we allocate in.
          if (allocated(topoDistance)) deallocate(topoDistance)
          if (allocated(toposForAtom)) deallocate(toposForAtom)
          if (allocated(nTopoPairsForAtom)) deallocate(nTopoPairsForAtom)
          if (allocated(skipsForLocalAtom)) deallocate(skipsForLocalAtom)
          if (allocated(nSkipsForLocalAtom)) deallocate(nSkipsForLocalAtom)
          if (allocated(skipsForRowAtom)) deallocate(skipsForRowAtom)
          if (allocated(nSkipsForRowAtom)) deallocate(nSkipsForRowAtom)
          if (allocated(mfactLocal)) deallocate(mfactLocal)
          if (allocated(mfactCol)) deallocate(mfactCol)
          if (allocated(mfactRow)) deallocate(mfactRow)
          if (allocated(groupListCol)) deallocate(groupListCol)     
          if (allocated(groupListRow)) deallocate(groupListRow)    
          if (allocated(groupStartCol)) deallocate(groupStartCol)
          if (allocated(groupStartRow)) deallocate(groupStartRow)     
          if (allocated(molMembershipList)) deallocate(molMembershipList)    
          if (allocated(excludes)) deallocate(excludes)

        end subroutine FreeSimGlobals

        pure function getNlocal() result(n)
          integer :: n
          n = nLocal
        end function getNlocal

        subroutine setHmatDangerousRcutValue(dangerWillRobinson)
          real(kind=dp), intent(in) :: dangerWillRobinson
          DangerRcut = dangerWillRobinson

          call checkBox()

          return
        end subroutine setHmatDangerousRcutValue
        
      end module simulation
      
