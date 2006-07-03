!!
!! Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
!!
!! The University of Notre Dame grants you ("Licensee") a
!! non-exclusive, royalty free, license to use, modify and
!! redistribute this software in source and binary code form, provided
!! that the following conditions are met:
!!
!! 1. Acknowledgement of the program authors must be made in any
!!    publication of scientific results based in part on use of the
!!    program.  An acceptable form of acknowledgement is citation of
!!    the article in which the program was described (Matthew
!!    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
!!    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
!!    Parallel Simulation Engine for Molecular Dynamics,"
!!    J. Comput. Chem. 26, pp. 252-271 (2005))
!!
!! 2. Redistributions of source code must retain the above copyright
!!    notice, this list of conditions and the following disclaimer.
!!
!! 3. Redistributions in binary form must reproduce the above copyright
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
#include "UseTheForce/DarkSide/fElectrostaticSummationMethod.h"

  type (simtype), public, save :: thisSim

  logical, save :: simulation_setup_complete = .false.

  integer, public, save :: nLocal, nGlobal
  integer, public, save :: nGroups, nGroupGlobal
  integer, public, save :: nExcludes_Global = 0
  integer, public, save :: nExcludes_Local = 0
  integer, allocatable, dimension(:,:), public :: excludesLocal
  integer, allocatable, dimension(:),   public :: excludesGlobal
  integer, allocatable, dimension(:),   public :: molMembershipList
  integer, allocatable, dimension(:),   public :: groupListRow
  integer, allocatable, dimension(:),   public :: groupStartRow
  integer, allocatable, dimension(:),   public :: groupListCol
  integer, allocatable, dimension(:),   public :: groupStartCol
  integer, allocatable, dimension(:),   public :: groupListLocal
  integer, allocatable, dimension(:),   public :: groupStartLocal
  integer, allocatable, dimension(:),   public :: nSkipsForAtom
  integer, allocatable, dimension(:,:), public :: skipsForAtom
  real(kind=dp), allocatable, dimension(:), public :: mfactRow
  real(kind=dp), allocatable, dimension(:), public :: mfactCol
  real(kind=dp), allocatable, dimension(:), public :: mfactLocal

  logical, allocatable, dimension(:) :: simHasAtypeMap
#ifdef IS_MPI
  logical, allocatable, dimension(:) :: simHasAtypeMapTemp
#endif

  real(kind=dp), public, dimension(3,3), save :: Hmat, HmatInv
  real(kind=dp), save :: DangerRcut
  logical, public, save :: boxIsOrthorhombic

  public :: SimulationSetup
  public :: getNlocal
  public :: setBox
  public :: checkBox
  public :: getDielect
  public :: SimUsesPBC

  public :: SimUsesDirectionalAtoms
  public :: SimUsesLennardJones
  public :: SimUsesElectrostatics
  public :: SimUsesCharges
  public :: SimUsesDipoles
  public :: SimUsesSticky
  public :: SimUsesStickyPower
  public :: SimUsesGayBerne
  public :: SimUsesEAM
  public :: SimUsesShapes
  public :: SimUsesFLARB
  public :: SimUsesRF
  public :: SimUsesSF
  public :: SimUsesSP
  public :: SimUsesBoxDipole
  public :: SimRequiresPrepairCalc
  public :: SimRequiresPostpairCalc
  public :: SimHasAtype
  public :: SimUsesSC
  public :: SimUsesMEAM
  public :: setHmatDangerousRcutValue

contains

  subroutine SimulationSetup(setThisSim, CnGlobal, CnLocal, c_idents, &
       CnLocalExcludes, CexcludesLocal, CnGlobalExcludes, CexcludesGlobal, &
       CmolMembership, Cmfact, CnGroups, CglobalGroupMembership, &
       status)    

    type (simtype) :: setThisSim
    integer, intent(inout) :: CnGlobal, CnLocal
    integer, dimension(CnLocal),intent(inout) :: c_idents

    integer :: CnLocalExcludes
    integer, dimension(2,CnLocalExcludes), intent(in) :: CexcludesLocal
    integer :: CnGlobalExcludes
    integer, dimension(CnGlobalExcludes), intent(in) :: CexcludesGlobal
    integer, dimension(CnGlobal),intent(in) :: CmolMembership 
    !!  Result status, success = 0, status = -1
    integer, intent(out) :: status
    integer :: i, j, me, thisStat, alloc_stat, myNode, id1, id2
    integer :: ia

    !! mass factors used for molecular cutoffs
    real ( kind = dp ), dimension(CnLocal) :: Cmfact
    integer, intent(in):: CnGroups
    integer, dimension(CnGlobal), intent(in):: CglobalGroupMembership
    integer :: maxSkipsForAtom, glPointer

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

    nExcludes_Global = CnGlobalExcludes
    nExcludes_Local = CnLocalExcludes

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

    allocate(c_idents_Row(nAtomsInRow),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1
       return
    endif

    allocate(c_idents_Col(nAtomsInCol),stat=alloc_stat)
    if (alloc_stat /= 0 ) then
       status = -1
       return
    endif

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

    !! free temporary ident arrays
    if (allocated(c_idents_Col)) then
       deallocate(c_idents_Col)
    end if
    if (allocated(c_idents_Row)) then
       deallocate(c_idents_Row)
    endif

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

    enddo

    do i = 1, nExcludes_Local
       excludesLocal(1,i) = CexcludesLocal(1,i)
       excludesLocal(2,i) = CexcludesLocal(2,i)
    enddo

#ifdef IS_MPI
    allocate(nSkipsForAtom(nAtomsInRow), stat=alloc_stat)
#else
    allocate(nSkipsForAtom(nLocal), stat=alloc_stat)
#endif
    if (alloc_stat /= 0 ) then
       thisStat = -1
       write(*,*) 'Could not allocate nSkipsForAtom array'
       return
    endif

    maxSkipsForAtom = 0
#ifdef IS_MPI
    do j = 1, nAtomsInRow
#else
       do j = 1, nLocal
#endif
          nSkipsForAtom(j) = 0
#ifdef IS_MPI
          id1 = AtomRowToGlobal(j)
#else 
          id1 = j
#endif
          do i = 1, nExcludes_Local
             if (excludesLocal(1,i) .eq. id1 ) then
                nSkipsForAtom(j) = nSkipsForAtom(j) + 1

                if (nSkipsForAtom(j) .gt. maxSkipsForAtom) then
                   maxSkipsForAtom = nSkipsForAtom(j)
                endif
             endif
             if (excludesLocal(2,i) .eq. id1 ) then
                nSkipsForAtom(j) = nSkipsForAtom(j) + 1

                if (nSkipsForAtom(j) .gt. maxSkipsForAtom) then
                   maxSkipsForAtom = nSkipsForAtom(j)
                endif
             endif
          end do
       enddo

#ifdef IS_MPI
       allocate(skipsForAtom(nAtomsInRow, maxSkipsForAtom), stat=alloc_stat)
#else
       allocate(skipsForAtom(nLocal, maxSkipsForAtom), stat=alloc_stat)
#endif
       if (alloc_stat /= 0 ) then
          write(*,*) 'Could not allocate skipsForAtom array'
          return
       endif

#ifdef IS_MPI
       do j = 1, nAtomsInRow
#else
          do j = 1, nLocal
#endif
             nSkipsForAtom(j) = 0
#ifdef IS_MPI
             id1 = AtomRowToGlobal(j)
#else 
             id1 = j
#endif
             do i = 1, nExcludes_Local
                if (excludesLocal(1,i) .eq. id1 ) then
                   nSkipsForAtom(j) = nSkipsForAtom(j) + 1
                   ! exclude lists have global ID's so this line is 
                   ! the same in MPI and non-MPI
                   id2 = excludesLocal(2,i)
                   skipsForAtom(j, nSkipsForAtom(j)) = id2
                endif
                if (excludesLocal(2, i) .eq. id1 ) then
                   nSkipsForAtom(j) = nSkipsForAtom(j) + 1
                   ! exclude lists have global ID's so this line is 
                   ! the same in MPI and non-MPI
                   id2 = excludesLocal(1,i)
                   skipsForAtom(j, nSkipsForAtom(j)) = id2
                endif
             end do
          enddo

          do i = 1, nExcludes_Global
             excludesGlobal(i) = CexcludesGlobal(i)
          enddo

          do i = 1, nGlobal
             molMemberShipList(i) = CmolMembership(i)
          enddo

         call createSimHasAtype(alloc_stat)
         if (alloc_stat /= 0) then
            status = -1
         end if
         
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

        function getDielect() result(dielect)
          real( kind = dp ) :: dielect
          dielect = thisSim%dielect
        end function getDielect

        function SimUsesPBC() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_PBC
        end function SimUsesPBC

        function SimUsesDirectionalAtoms() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_dipoles .or. thisSim%SIM_uses_Sticky .or. &
               thisSim%SIM_uses_StickyPower .or. &
               thisSim%SIM_uses_GayBerne .or. thisSim%SIM_uses_Shapes
        end function SimUsesDirectionalAtoms

        function SimUsesLennardJones() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_LennardJones
        end function SimUsesLennardJones

        function SimUsesElectrostatics() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_Electrostatics
        end function SimUsesElectrostatics

        function SimUsesCharges() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_Charges
        end function SimUsesCharges

        function SimUsesDipoles() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_Dipoles
        end function SimUsesDipoles

        function SimUsesSticky() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_Sticky
        end function SimUsesSticky

        function SimUsesStickyPower() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_StickyPower
        end function SimUsesStickyPower

        function SimUsesGayBerne() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_GayBerne
        end function SimUsesGayBerne

        function SimUsesEAM() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_EAM
        end function SimUsesEAM


        function SimUsesSC() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_SC
        end function SimUsesSC

        function SimUsesMEAM() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_MEAM
        end function SimUsesMEAM


        function SimUsesShapes() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_Shapes
        end function SimUsesShapes

        function SimUsesFLARB() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_FLARB
        end function SimUsesFLARB

        function SimUsesRF() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_RF
        end function SimUsesRF

        function SimUsesSF() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_SF
        end function SimUsesSF

        function SimUsesSP() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_SP
        end function SimUsesSP

        function SimUsesBoxDipole() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_BoxDipole
        end function SimUsesBoxDipole

        function SimRequiresPrepairCalc() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_EAM .or. thisSim%SIM_uses_SC &
               .or. thisSim%SIM_uses_MEAM
        end function SimRequiresPrepairCalc
        
        function SimRequiresPostpairCalc() result(doesit)
          logical :: doesit
          doesit = thisSim%SIM_uses_RF .or. thisSim%SIM_uses_SF &
               .or. thisSim%SIM_uses_SP .or. thisSim%SIM_uses_BoxDipole
        end function SimRequiresPostpairCalc

        ! Function returns true if the simulation has this atype
        function SimHasAtype(thisAtype) result(doesit)
          logical :: doesit
          integer :: thisAtype
          doesit = .false.
          if(.not.allocated(SimHasAtypeMap)) return

          doesit = SimHasAtypeMap(thisAtype)
             
        end function SimHasAtype

        subroutine createSimHasAtype(status)
          integer, intent(out) :: status
          integer :: alloc_stat
          integer :: me_i
          integer :: mpiErrors
          integer :: nAtypes
          status = 0

          nAtypes = getSize(atypes)
          ! Setup logical map for atypes in simulation
          if (.not.allocated(SimHasAtypeMap)) then
             allocate(SimHasAtypeMap(nAtypes),stat=alloc_stat)
             if (alloc_stat /= 0 ) then
                status = -1
                return
             end if
             SimHasAtypeMap = .false.
          end if
          
          ! Loop through the local atoms and grab the atypes present         
          do me_i = 1,nLocal
             SimHasAtypeMap(atid(me_i)) = .true.
          end do
          ! For MPI, we need to know all possible atypes present in 
          ! simulation on all processors. Use LOR operation to set map.
#ifdef IS_MPI
          if (.not.allocated(SimHasAtypeMapTemp)) then
             allocate(SimHasAtypeMapTemp(nAtypes),stat=alloc_stat)
             if (alloc_stat /= 0 ) then
                status = -1
                return
             end if
          end if
          call mpi_allreduce(SimHasAtypeMap, SimHasAtypeMaptemp, nAtypes, &
               mpi_logical, MPI_LOR, mpi_comm_world, mpiErrors)
          simHasAtypeMap = simHasAtypeMapTemp
          deallocate(simHasAtypeMapTemp)
#endif          
        end subroutine createSimHasAtype
        
       subroutine InitializeSimGlobals(thisStat)
          integer, intent(out) :: thisStat
          integer :: alloc_stat

          thisStat = 0

          call FreeSimGlobals()    

          allocate(excludesLocal(2,nExcludes_Local), stat=alloc_stat)
          if (alloc_stat /= 0 ) then
             thisStat = -1
             return
          endif

          allocate(excludesGlobal(nExcludes_Global), stat=alloc_stat)
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

          if (allocated(skipsForAtom)) deallocate(skipsForAtom)
          if (allocated(nSkipsForAtom)) deallocate(nSkipsForAtom)
          if (allocated(mfactLocal)) deallocate(mfactLocal)
          if (allocated(mfactCol)) deallocate(mfactCol)
          if (allocated(mfactRow)) deallocate(mfactRow)
          if (allocated(groupListCol)) deallocate(groupListCol)     
          if (allocated(groupListRow)) deallocate(groupListRow)    
          if (allocated(groupStartCol)) deallocate(groupStartCol)
          if (allocated(groupStartRow)) deallocate(groupStartRow)     
          if (allocated(molMembershipList)) deallocate(molMembershipList)     
          if (allocated(excludesGlobal)) deallocate(excludesGlobal)
          if (allocated(excludesLocal)) deallocate(excludesLocal)

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
