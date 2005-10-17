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

!! doForces.F90
!! module doForces
!! Calculates Long Range forces.

!! @author Charles F. Vardeman II
!! @author Matthew Meineke
!! @version $Id: doForces.F90,v 1.59 2005-10-17 19:12:34 gezelter Exp $, $Date: 2005-10-17 19:12:34 $, $Name: not supported by cvs2svn $, $Revision: 1.59 $


module doForces
  use force_globals
  use simulation
  use definitions
  use atype_module
  use switcheroo
  use neighborLists  
  use lj
  use sticky
  use electrostatic_module
  use reaction_field_module
  use gayberne
  use shapes
  use vector_class
  use eam
  use status
#ifdef IS_MPI
  use mpiSimulation
#endif

  implicit none
  PRIVATE

#define __FORTRAN90
#include "UseTheForce/fSwitchingFunction.h"
#include "UseTheForce/fCutoffPolicy.h"
#include "UseTheForce/DarkSide/fInteractionMap.h"
#include "UseTheForce/DarkSide/fElectrostaticSummationMethod.h"


  INTEGER, PARAMETER:: PREPAIR_LOOP = 1
  INTEGER, PARAMETER:: PAIR_LOOP    = 2

  logical, save :: haveNeighborList = .false.
  logical, save :: haveSIMvariables = .false.
  logical, save :: haveSaneForceField = .false.
  logical, save :: haveInteractionHash = .false.
  logical, save :: haveGtypeCutoffMap = .false.
  logical, save :: haveDefaultCutoffs = .false.
  logical, save :: haveRlist = .false.

  logical, save :: FF_uses_DirectionalAtoms
  logical, save :: FF_uses_Dipoles
  logical, save :: FF_uses_GayBerne
  logical, save :: FF_uses_EAM

  logical, save :: SIM_uses_DirectionalAtoms
  logical, save :: SIM_uses_EAM
  logical, save :: SIM_requires_postpair_calc
  logical, save :: SIM_requires_prepair_calc
  logical, save :: SIM_uses_PBC

  integer, save :: electrostaticSummationMethod

  public :: init_FF
  public :: setDefaultCutoffs
  public :: do_force_loop
  public :: createInteractionHash
  public :: createGtypeCutoffMap
  public :: getStickyCut
  public :: getStickyPowerCut
  public :: getGayBerneCut
  public :: getEAMCut
  public :: getShapeCut

#ifdef PROFILE
  public :: getforcetime
  real, save :: forceTime = 0
  real :: forceTimeInitial, forceTimeFinal
  integer :: nLoops
#endif
  
  !! Variables for cutoff mapping and interaction mapping
  ! Bit hash to determine pair-pair interactions.
  integer, dimension(:,:), allocatable :: InteractionHash
  real(kind=dp), dimension(:), allocatable :: atypeMaxCutoff
  real(kind=dp), dimension(:), allocatable, target :: groupMaxCutoffRow
  real(kind=dp), dimension(:), pointer :: groupMaxCutoffCol

  integer, dimension(:), allocatable, target :: groupToGtypeRow
  integer, dimension(:), pointer :: groupToGtypeCol => null()

  real(kind=dp), dimension(:), allocatable,target :: gtypeMaxCutoffRow
  real(kind=dp), dimension(:), pointer :: gtypeMaxCutoffCol
  type ::gtypeCutoffs
     real(kind=dp) :: rcut 
     real(kind=dp) :: rcutsq 
     real(kind=dp) :: rlistsq
  end type gtypeCutoffs
  type(gtypeCutoffs), dimension(:,:), allocatable :: gtypeCutoffMap

  integer, save :: cutoffPolicy = TRADITIONAL_CUTOFF_POLICY
  real(kind=dp),save :: defaultRcut, defaultRsw, defaultRlist
  real(kind=dp),save :: listSkin
  
contains

  subroutine createInteractionHash(status)
    integer :: nAtypes
    integer, intent(out) :: status
    integer :: i
    integer :: j
    integer :: iHash
    !! Test Types
    logical :: i_is_LJ
    logical :: i_is_Elect
    logical :: i_is_Sticky
    logical :: i_is_StickyP
    logical :: i_is_GB
    logical :: i_is_EAM
    logical :: i_is_Shape
    logical :: j_is_LJ
    logical :: j_is_Elect
    logical :: j_is_Sticky
    logical :: j_is_StickyP
    logical :: j_is_GB
    logical :: j_is_EAM
    logical :: j_is_Shape
    real(kind=dp) :: myRcut

    status = 0   

    if (.not. associated(atypes)) then
       call handleError("atype", "atypes was not present before call of createInteractionHash!")
       status = -1
       return
    endif
    
    nAtypes = getSize(atypes)
    
    if (nAtypes == 0) then
       status = -1
       return
    end if

    if (.not. allocated(InteractionHash)) then
       allocate(InteractionHash(nAtypes,nAtypes))
    else
       deallocate(InteractionHash)
       allocate(InteractionHash(nAtypes,nAtypes))
    endif

    if (.not. allocated(atypeMaxCutoff)) then
       allocate(atypeMaxCutoff(nAtypes))
    else
       deallocate(atypeMaxCutoff)
       allocate(atypeMaxCutoff(nAtypes))
    endif
        
    do i = 1, nAtypes
       call getElementProperty(atypes, i, "is_LennardJones", i_is_LJ)
       call getElementProperty(atypes, i, "is_Electrostatic", i_is_Elect)
       call getElementProperty(atypes, i, "is_Sticky", i_is_Sticky)
       call getElementProperty(atypes, i, "is_StickyPower", i_is_StickyP)
       call getElementProperty(atypes, i, "is_GayBerne", i_is_GB)
       call getElementProperty(atypes, i, "is_EAM", i_is_EAM)
       call getElementProperty(atypes, i, "is_Shape", i_is_Shape)

       do j = i, nAtypes

          iHash = 0
          myRcut = 0.0_dp

          call getElementProperty(atypes, j, "is_LennardJones", j_is_LJ)
          call getElementProperty(atypes, j, "is_Electrostatic", j_is_Elect)
          call getElementProperty(atypes, j, "is_Sticky", j_is_Sticky)
          call getElementProperty(atypes, j, "is_StickyPower", j_is_StickyP)
          call getElementProperty(atypes, j, "is_GayBerne", j_is_GB)
          call getElementProperty(atypes, j, "is_EAM", j_is_EAM)
          call getElementProperty(atypes, j, "is_Shape", j_is_Shape)

          if (i_is_LJ .and. j_is_LJ) then
             iHash = ior(iHash, LJ_PAIR)            
          endif
          
          if (i_is_Elect .and. j_is_Elect) then
             iHash = ior(iHash, ELECTROSTATIC_PAIR)
          endif
          
          if (i_is_Sticky .and. j_is_Sticky) then
             iHash = ior(iHash, STICKY_PAIR)
          endif

          if (i_is_StickyP .and. j_is_StickyP) then
             iHash = ior(iHash, STICKYPOWER_PAIR)
          endif

          if (i_is_EAM .and. j_is_EAM) then
             iHash = ior(iHash, EAM_PAIR)
          endif

          if (i_is_GB .and. j_is_GB) iHash = ior(iHash, GAYBERNE_PAIR)
          if (i_is_GB .and. j_is_LJ) iHash = ior(iHash, GAYBERNE_LJ)
          if (i_is_LJ .and. j_is_GB) iHash = ior(iHash, GAYBERNE_LJ)

          if (i_is_Shape .and. j_is_Shape) iHash = ior(iHash, SHAPE_PAIR)
          if (i_is_Shape .and. j_is_LJ) iHash = ior(iHash, SHAPE_LJ)
          if (i_is_LJ .and. j_is_Shape) iHash = ior(iHash, SHAPE_LJ)


          InteractionHash(i,j) = iHash
          InteractionHash(j,i) = iHash

       end do

    end do

    haveInteractionHash = .true.
  end subroutine createInteractionHash

  subroutine createGtypeCutoffMap(stat)

    integer, intent(out), optional :: stat
    logical :: i_is_LJ
    logical :: i_is_Elect
    logical :: i_is_Sticky
    logical :: i_is_StickyP
    logical :: i_is_GB
    logical :: i_is_EAM
    logical :: i_is_Shape
    logical :: GtypeFound

    integer :: myStatus, nAtypes,  i, j, istart, iend, jstart, jend
    integer :: n_in_i, me_i, ia, g, atom1, ja, n_in_j,me_j
    integer :: nGroupsInRow
    integer :: nGroupsInCol
    integer :: nGroupTypesRow,nGroupTypesCol
    real(kind=dp):: thisSigma, bigSigma, thisRcut, tradRcut, tol, skin
    real(kind=dp) :: biggestAtypeCutoff

    stat = 0
    if (.not. haveInteractionHash) then 
       call createInteractionHash(myStatus)      
       if (myStatus .ne. 0) then
          write(default_error, *) 'createInteractionHash failed in doForces!'
          stat = -1
          return
       endif
    endif
#ifdef IS_MPI
    nGroupsInRow = getNgroupsInRow(plan_group_row)
    nGroupsInCol = getNgroupsInCol(plan_group_col)
#endif
    nAtypes = getSize(atypes)
! Set all of the initial cutoffs to zero.
    atypeMaxCutoff = 0.0_dp
    do i = 1, nAtypes
       if (SimHasAtype(i)) then     
          call getElementProperty(atypes, i, "is_LennardJones", i_is_LJ)
          call getElementProperty(atypes, i, "is_Electrostatic", i_is_Elect)
          call getElementProperty(atypes, i, "is_Sticky", i_is_Sticky)
          call getElementProperty(atypes, i, "is_StickyPower", i_is_StickyP)
          call getElementProperty(atypes, i, "is_GayBerne", i_is_GB)
          call getElementProperty(atypes, i, "is_EAM", i_is_EAM)
          call getElementProperty(atypes, i, "is_Shape", i_is_Shape)
          
 
          if (haveDefaultCutoffs) then
             atypeMaxCutoff(i) = defaultRcut
          else
             if (i_is_LJ) then           
                thisRcut = getSigma(i) * 2.5_dp
                if (thisRCut .gt. atypeMaxCutoff(i)) atypeMaxCutoff(i) = thisRCut
             endif
             if (i_is_Elect) then
                thisRcut = defaultRcut
                if (thisRCut .gt. atypeMaxCutoff(i)) atypeMaxCutoff(i) = thisRCut
             endif
             if (i_is_Sticky) then
                thisRcut = getStickyCut(i)
                if (thisRCut .gt. atypeMaxCutoff(i)) atypeMaxCutoff(i) = thisRCut
             endif
             if (i_is_StickyP) then
                thisRcut = getStickyPowerCut(i)
                if (thisRCut .gt. atypeMaxCutoff(i)) atypeMaxCutoff(i) = thisRCut
             endif
             if (i_is_GB) then
                thisRcut = getGayBerneCut(i)
                if (thisRCut .gt. atypeMaxCutoff(i)) atypeMaxCutoff(i) = thisRCut
             endif
             if (i_is_EAM) then
                thisRcut = getEAMCut(i)
                if (thisRCut .gt. atypeMaxCutoff(i)) atypeMaxCutoff(i) = thisRCut
             endif
             if (i_is_Shape) then
                thisRcut = getShapeCut(i)
                if (thisRCut .gt. atypeMaxCutoff(i)) atypeMaxCutoff(i) = thisRCut
             endif
          endif
          
          
          if (atypeMaxCutoff(i).gt.biggestAtypeCutoff) then
             biggestAtypeCutoff = atypeMaxCutoff(i)
          endif

       endif
    enddo
  

    
    istart = 1
    jstart = 1
#ifdef IS_MPI
    iend = nGroupsInRow
    jend = nGroupsInCol
#else
    iend = nGroups 
    jend = nGroups
#endif
    
    !! allocate the groupToGtype and gtypeMaxCutoff here.
    if(.not.allocated(groupToGtypeRow)) then
     !  allocate(groupToGtype(iend))
       allocate(groupToGtypeRow(iend))
    else
       deallocate(groupToGtypeRow)
       allocate(groupToGtypeRow(iend))
    endif
    if(.not.allocated(groupMaxCutoffRow)) then
       allocate(groupMaxCutoffRow(iend))
    else
       deallocate(groupMaxCutoffRow)
       allocate(groupMaxCutoffRow(iend))
    end if

    if(.not.allocated(gtypeMaxCutoffRow)) then
       allocate(gtypeMaxCutoffRow(iend))
    else
       deallocate(gtypeMaxCutoffRow)
       allocate(gtypeMaxCutoffRow(iend))
    endif


#ifdef IS_MPI
       ! We only allocate new storage if we are in MPI because Ncol /= Nrow
    if(.not.associated(groupToGtypeCol)) then
       allocate(groupToGtypeCol(jend))
    else
       deallocate(groupToGtypeCol)
       allocate(groupToGtypeCol(jend))
    end if

    if(.not.associated(groupToGtypeCol)) then
       allocate(groupToGtypeCol(jend))
    else
       deallocate(groupToGtypeCol)
       allocate(groupToGtypeCol(jend))
    end if
    if(.not.associated(gtypeMaxCutoffCol)) then
       allocate(gtypeMaxCutoffCol(jend))
    else
       deallocate(gtypeMaxCutoffCol)       
       allocate(gtypeMaxCutoffCol(jend))
    end if

       groupMaxCutoffCol = 0.0_dp
       gtypeMaxCutoffCol = 0.0_dp

#endif
       groupMaxCutoffRow = 0.0_dp
       gtypeMaxCutoffRow = 0.0_dp


    !! first we do a single loop over the cutoff groups to find the
    !! largest cutoff for any atypes present in this group.  We also
    !! create gtypes at this point.
    
    tol = 1.0d-6
    nGroupTypesRow = 0

    do i = istart, iend       
       n_in_i = groupStartRow(i+1) - groupStartRow(i)
       groupMaxCutoffRow(i) = 0.0_dp
       do ia = groupStartRow(i), groupStartRow(i+1)-1
          atom1 = groupListRow(ia)
#ifdef IS_MPI
          me_i = atid_row(atom1)
#else
          me_i = atid(atom1)
#endif          
          if (atypeMaxCutoff(me_i).gt.groupMaxCutoffRow(i)) then 
             groupMaxCutoffRow(i)=atypeMaxCutoff(me_i)
          endif          
       enddo

       if (nGroupTypesRow.eq.0) then
          nGroupTypesRow = nGroupTypesRow + 1
          gtypeMaxCutoffRow(nGroupTypesRow) = groupMaxCutoffRow(i)
          groupToGtypeRow(i) = nGroupTypesRow
       else
          GtypeFound = .false.
          do g = 1, nGroupTypesRow
             if ( abs(groupMaxCutoffRow(i) - gtypeMaxCutoffRow(g)).lt.tol) then
                groupToGtypeRow(i) = g
                GtypeFound = .true.
             endif
          enddo
          if (.not.GtypeFound) then             
             nGroupTypesRow = nGroupTypesRow + 1
             gtypeMaxCutoffRow(nGroupTypesRow) = groupMaxCutoffRow(i)
             groupToGtypeRow(i) = nGroupTypesRow
          endif
       endif
    enddo    

#ifdef IS_MPI
    do j = jstart, jend       
       n_in_j = groupStartCol(j+1) - groupStartCol(j)
       groupMaxCutoffCol(j) = 0.0_dp
       do ja = groupStartCol(j), groupStartCol(j+1)-1
          atom1 = groupListCol(ja)

          me_j = atid_col(atom1)

          if (atypeMaxCutoff(me_j).gt.groupMaxCutoffCol(j)) then 
             groupMaxCutoffCol(j)=atypeMaxCutoff(me_j)
          endif          
       enddo

       if (nGroupTypesCol.eq.0) then
          nGroupTypesCol = nGroupTypesCol + 1
          gtypeMaxCutoffCol(nGroupTypesCol) = groupMaxCutoffCol(j)
          groupToGtypeCol(j) = nGroupTypesCol
       else
          GtypeFound = .false.
          do g = 1, nGroupTypesCol
             if ( abs(groupMaxCutoffCol(j) - gtypeMaxCutoffCol(g)).lt.tol) then
                groupToGtypeCol(j) = g
                GtypeFound = .true.
             endif
          enddo
          if (.not.GtypeFound) then             
             nGroupTypesCol = nGroupTypesCol + 1
             gtypeMaxCutoffCol(nGroupTypesCol) = groupMaxCutoffCol(j)
             groupToGtypeCol(j) = nGroupTypesCol
          endif
       endif
    enddo    

#else
! Set pointers to information we just found
    nGroupTypesCol = nGroupTypesRow
    groupToGtypeCol => groupToGtypeRow
    gtypeMaxCutoffCol => gtypeMaxCutoffRow
    groupMaxCutoffCol => groupMaxCutoffRow
#endif





    !! allocate the gtypeCutoffMap here.
    allocate(gtypeCutoffMap(nGroupTypesRow,nGroupTypesCol))
    !! then we do a double loop over all the group TYPES to find the cutoff
    !! map between groups of two types
    tradRcut = max(maxval(gtypeMaxCutoffRow),maxval(gtypeMaxCutoffCol))

    do i = 1, nGroupTypesRow
       do j = 1, nGroupTypesCol
       
          select case(cutoffPolicy)
          case(TRADITIONAL_CUTOFF_POLICY)
             thisRcut = tradRcut
          case(MIX_CUTOFF_POLICY)
             thisRcut = 0.5_dp * (gtypeMaxCutoffRow(i) + gtypeMaxCutoffCol(j))
          case(MAX_CUTOFF_POLICY)
             thisRcut = max(gtypeMaxCutoffRow(i), gtypeMaxCutoffCol(j))
          case default
             call handleError("createGtypeCutoffMap", "Unknown Cutoff Policy")
             return
          end select
          gtypeCutoffMap(i,j)%rcut = thisRcut
          gtypeCutoffMap(i,j)%rcutsq = thisRcut*thisRcut
          skin = defaultRlist - defaultRcut
          listSkin = skin ! set neighbor list skin thickness
          gtypeCutoffMap(i,j)%rlistsq = (thisRcut + skin)**2

          ! sanity check

          if (haveDefaultCutoffs) then
             if (abs(gtypeCutoffMap(i,j)%rcut - defaultRcut).gt.0.0001) then
                call handleError("createGtypeCutoffMap", "user-specified rCut does not match computed group Cutoff")
             endif
          endif
       enddo
    enddo
    if(allocated(gtypeMaxCutoffRow)) deallocate(gtypeMaxCutoffRow)
    if(allocated(groupMaxCutoffRow)) deallocate(groupMaxCutoffRow)
    if(allocated(atypeMaxCutoff)) deallocate(atypeMaxCutoff)
#ifdef IS_MPI
    if(associated(groupMaxCutoffCol)) deallocate(groupMaxCutoffCol)
    if(associated(gtypeMaxCutoffCol)) deallocate(gtypeMaxCutoffCol)
#endif
    groupMaxCutoffCol => null()
    gtypeMaxCutoffCol => null()
    
    haveGtypeCutoffMap = .true.
   end subroutine createGtypeCutoffMap

   subroutine setDefaultCutoffs(defRcut, defRsw, defRlist, cutPolicy)
     real(kind=dp),intent(in) :: defRcut, defRsw, defRlist
     integer, intent(in) :: cutPolicy

     defaultRcut = defRcut
     defaultRsw = defRsw
     defaultRlist = defRlist
     cutoffPolicy = cutPolicy

     haveDefaultCutoffs = .true.
   end subroutine setDefaultCutoffs

   subroutine setCutoffPolicy(cutPolicy)

     integer, intent(in) :: cutPolicy
     cutoffPolicy = cutPolicy
     call createGtypeCutoffMap()
   end subroutine setCutoffPolicy
    
     
  subroutine setSimVariables()
    SIM_uses_DirectionalAtoms = SimUsesDirectionalAtoms()
    SIM_uses_EAM = SimUsesEAM()
    SIM_requires_postpair_calc = SimRequiresPostpairCalc()
    SIM_requires_prepair_calc = SimRequiresPrepairCalc()
    SIM_uses_PBC = SimUsesPBC()

    haveSIMvariables = .true.

    return
  end subroutine setSimVariables

  subroutine doReadyCheck(error)
    integer, intent(out) :: error

    integer :: myStatus

    error = 0

    if (.not. haveInteractionHash) then       
       myStatus = 0       
       call createInteractionHash(myStatus)       
       if (myStatus .ne. 0) then
          write(default_error, *) 'createInteractionHash failed in doForces!'
          error = -1
          return
       endif
    endif

    if (.not. haveGtypeCutoffMap) then        
       myStatus = 0       
       call createGtypeCutoffMap(myStatus)       
       if (myStatus .ne. 0) then
          write(default_error, *) 'createGtypeCutoffMap failed in doForces!'
          error = -1
          return
       endif
    endif

    if (.not. haveSIMvariables) then
       call setSimVariables()
    endif

  !  if (.not. haveRlist) then
  !     write(default_error, *) 'rList has not been set in doForces!'
  !     error = -1
  !     return
  !  endif

    if (.not. haveNeighborList) then
       write(default_error, *) 'neighbor list has not been initialized in doForces!'
       error = -1
       return
    end if

    if (.not. haveSaneForceField) then
       write(default_error, *) 'Force Field is not sane in doForces!'
       error = -1
       return
    end if

#ifdef IS_MPI
    if (.not. isMPISimSet()) then
       write(default_error,*) "ERROR: mpiSimulation has not been initialized!"
       error = -1
       return
    endif
#endif
    return
  end subroutine doReadyCheck


  subroutine init_FF(thisESM, thisStat)

    integer, intent(in) :: thisESM
    integer, intent(out) :: thisStat   
    integer :: my_status, nMatches
    integer, pointer :: MatchList(:) => null()
    real(kind=dp) :: rcut, rrf, rt, dielect

    !! assume things are copacetic, unless they aren't
    thisStat = 0

    electrostaticSummationMethod = thisESM

    !! init_FF is called *after* all of the atom types have been 
    !! defined in atype_module using the new_atype subroutine.
    !!
    !! this will scan through the known atypes and figure out what
    !! interactions are used by the force field.     

    FF_uses_DirectionalAtoms = .false.
    FF_uses_Dipoles = .false.
    FF_uses_GayBerne = .false.
    FF_uses_EAM = .false.

    call getMatchingElementList(atypes, "is_Directional", .true., &
         nMatches, MatchList)
    if (nMatches .gt. 0) FF_uses_DirectionalAtoms = .true.

    call getMatchingElementList(atypes, "is_Dipole", .true., &
         nMatches, MatchList)
    if (nMatches .gt. 0) FF_uses_Dipoles = .true.
    
    call getMatchingElementList(atypes, "is_GayBerne", .true., &
         nMatches, MatchList)
    if (nMatches .gt. 0) FF_uses_GayBerne = .true.

    call getMatchingElementList(atypes, "is_EAM", .true., nMatches, MatchList)
    if (nMatches .gt. 0) FF_uses_EAM = .true.


    haveSaneForceField = .true.

    !! check to make sure the reaction field setting makes sense

    if (FF_uses_Dipoles) then
       if (electrostaticSummationMethod == REACTION_FIELD) then
          dielect = getDielect()
          call initialize_rf(dielect)
       endif
    else
       if (electrostaticSummationMethod == REACTION_FIELD) then
          write(default_error,*) 'Using Reaction Field with no dipoles?  Huh?'
          thisStat = -1
          haveSaneForceField = .false.
          return
       endif
    endif

    if (FF_uses_EAM) then
       call init_EAM_FF(my_status) 
       if (my_status /= 0) then
          write(default_error, *) "init_EAM_FF returned a bad status"
          thisStat = -1
          haveSaneForceField = .false.
          return
       end if
    endif

    if (.not. haveNeighborList) then
       !! Create neighbor lists
       call expandNeighborList(nLocal, my_status)
       if (my_Status /= 0) then
          write(default_error,*) "SimSetup: ExpandNeighborList returned error."
          thisStat = -1
          return
       endif
       haveNeighborList = .true.
    endif

  end subroutine init_FF


  !! Does force loop over i,j pairs. Calls do_pair to calculates forces.
  !------------------------------------------------------------->
  subroutine do_force_loop(q, q_group, A, eFrame, f, t, tau, pot, &
       do_pot_c, do_stress_c, error)
    !! Position array provided by C, dimensioned by getNlocal
    real ( kind = dp ), dimension(3, nLocal) :: q
    !! molecular center-of-mass position array
    real ( kind = dp ), dimension(3, nGroups) :: q_group
    !! Rotation Matrix for each long range particle in simulation.
    real( kind = dp), dimension(9, nLocal) :: A    
    !! Unit vectors for dipoles (lab frame)
    real( kind = dp ), dimension(9,nLocal) :: eFrame
    !! Force array provided by C, dimensioned by getNlocal
    real ( kind = dp ), dimension(3,nLocal) :: f
    !! Torsion array provided by C, dimensioned by getNlocal
    real( kind = dp ), dimension(3,nLocal) :: t    

    !! Stress Tensor
    real( kind = dp), dimension(9) :: tau   
    real ( kind = dp ),dimension(LR_POT_TYPES) :: pot
    logical ( kind = 2) :: do_pot_c, do_stress_c
    logical :: do_pot
    logical :: do_stress
    logical :: in_switching_region
#ifdef IS_MPI 
    real( kind = DP ), dimension(LR_POT_TYPES) :: pot_local
    integer :: nAtomsInRow
    integer :: nAtomsInCol
    integer :: nprocs
    integer :: nGroupsInRow
    integer :: nGroupsInCol
#endif
    integer :: natoms    
    logical :: update_nlist   
    integer :: i, j, jstart, jend, jnab
    integer :: istart, iend
    integer :: ia, jb, atom1, atom2
    integer :: nlist
    real( kind = DP ) :: ratmsq, rgrpsq, rgrp, vpair, vij
    real( kind = DP ) :: sw, dswdr, swderiv, mf
    real(kind=dp),dimension(3) :: d_atm, d_grp, fpair, fij
    real(kind=dp) :: rfpot, mu_i, virial
    integer :: me_i, me_j, n_in_i, n_in_j
    logical :: is_dp_i
    integer :: neighborListSize
    integer :: listerror, error
    integer :: localError
    integer :: propPack_i, propPack_j
    integer :: loopStart, loopEnd, loop
    integer :: iHash
  

    !! initialize local variables  

#ifdef IS_MPI
    pot_local = 0.0_dp
    nAtomsInRow   = getNatomsInRow(plan_atom_row)
    nAtomsInCol   = getNatomsInCol(plan_atom_col)
    nGroupsInRow  = getNgroupsInRow(plan_group_row)
    nGroupsInCol  = getNgroupsInCol(plan_group_col)
#else
    natoms = nlocal
#endif

    call doReadyCheck(localError)
    if ( localError .ne. 0 ) then
       call handleError("do_force_loop", "Not Initialized")
       error = -1
       return
    end if
    call zero_work_arrays()

    do_pot = do_pot_c
    do_stress = do_stress_c

    ! Gather all information needed by all force loops:

#ifdef IS_MPI    

    call gather(q, q_Row, plan_atom_row_3d)
    call gather(q, q_Col, plan_atom_col_3d)

    call gather(q_group, q_group_Row, plan_group_row_3d)
    call gather(q_group, q_group_Col, plan_group_col_3d)

    if (FF_UsesDirectionalAtoms() .and. SIM_uses_DirectionalAtoms) then
       call gather(eFrame, eFrame_Row, plan_atom_row_rotation)
       call gather(eFrame, eFrame_Col, plan_atom_col_rotation)

       call gather(A, A_Row, plan_atom_row_rotation)
       call gather(A, A_Col, plan_atom_col_rotation)
    endif

#endif

    !! Begin force loop timing:
#ifdef PROFILE
    call cpu_time(forceTimeInitial)
    nloops = nloops + 1
#endif

    loopEnd = PAIR_LOOP
    if (FF_RequiresPrepairCalc() .and. SIM_requires_prepair_calc) then
       loopStart = PREPAIR_LOOP
    else
       loopStart = PAIR_LOOP
    endif

    do loop = loopStart, loopEnd

       ! See if we need to update neighbor lists
       ! (but only on the first time through):
       if (loop .eq. loopStart) then
#ifdef IS_MPI
          call checkNeighborList(nGroupsInRow, q_group_row, listSkin, &
               update_nlist)
#else
          call checkNeighborList(nGroups, q_group, listSkin, &
               update_nlist)
#endif
       endif

       if (update_nlist) then
          !! save current configuration and construct neighbor list
#ifdef IS_MPI
          call saveNeighborList(nGroupsInRow, q_group_row)
#else
          call saveNeighborList(nGroups, q_group)
#endif         
          neighborListSize = size(list)
          nlist = 0
       endif

       istart = 1
#ifdef IS_MPI
       iend = nGroupsInRow
#else
       iend = nGroups - 1
#endif
       outer: do i = istart, iend

          if (update_nlist) point(i) = nlist + 1

          n_in_i = groupStartRow(i+1) - groupStartRow(i)

          if (update_nlist) then
#ifdef IS_MPI
             jstart = 1
             jend = nGroupsInCol
#else
             jstart = i+1
             jend = nGroups
#endif
          else             
             jstart = point(i)
             jend = point(i+1) - 1
             ! make sure group i has neighbors
             if (jstart .gt. jend) cycle outer
          endif

          do jnab = jstart, jend
             if (update_nlist) then
                j = jnab
             else
                j = list(jnab)
             endif

#ifdef IS_MPI
             me_j = atid_col(j)
             call get_interatomic_vector(q_group_Row(:,i), &
                  q_group_Col(:,j), d_grp, rgrpsq)
#else
             me_j = atid(j)
             call get_interatomic_vector(q_group(:,i), &
                  q_group(:,j), d_grp, rgrpsq)
#endif      

             if (rgrpsq < gtypeCutoffMap(groupToGtypeRow(i),groupToGtypeCol(j))%rListsq) then
                if (update_nlist) then
                   nlist = nlist + 1

                   if (nlist > neighborListSize) then
#ifdef IS_MPI                 
                      call expandNeighborList(nGroupsInRow, listerror)
#else
                      call expandNeighborList(nGroups, listerror)
#endif
                      if (listerror /= 0) then
                         error = -1
                         write(DEFAULT_ERROR,*) "ERROR: nlist > list size and max allocations exceeded."
                         return
                      end if
                      neighborListSize = size(list)
                   endif

                   list(nlist) = j
                endif

                if (loop .eq. PAIR_LOOP) then
                   vij = 0.0d0
                   fij(1:3) = 0.0d0
                endif

                call get_switch(rgrpsq, sw, dswdr, rgrp, group_switch, &
                     in_switching_region)

                n_in_j = groupStartCol(j+1) - groupStartCol(j)

                do ia = groupStartRow(i), groupStartRow(i+1)-1

                   atom1 = groupListRow(ia)

                   inner: do jb = groupStartCol(j), groupStartCol(j+1)-1

                      atom2 = groupListCol(jb)

                      if (skipThisPair(atom1, atom2)) cycle inner

                      if ((n_in_i .eq. 1).and.(n_in_j .eq. 1)) then
                         d_atm(1:3) = d_grp(1:3)
                         ratmsq = rgrpsq
                      else
#ifdef IS_MPI
                         call get_interatomic_vector(q_Row(:,atom1), &
                              q_Col(:,atom2), d_atm, ratmsq)
#else
                         call get_interatomic_vector(q(:,atom1), &
                              q(:,atom2), d_atm, ratmsq)
#endif
                      endif

                      if (loop .eq. PREPAIR_LOOP) then
#ifdef IS_MPI                      
                         call do_prepair(atom1, atom2, ratmsq, d_atm, sw, &
                              rgrpsq, d_grp, do_pot, do_stress, &
                              eFrame, A, f, t, pot_local)
#else
                         call do_prepair(atom1, atom2, ratmsq, d_atm, sw, &
                              rgrpsq, d_grp, do_pot, do_stress, &
                              eFrame, A, f, t, pot)
#endif                                               
                      else
#ifdef IS_MPI                      
                         call do_pair(atom1, atom2, ratmsq, d_atm, sw, &
                              do_pot, &
                              eFrame, A, f, t, pot_local, vpair, fpair)
#else
                         call do_pair(atom1, atom2, ratmsq, d_atm, sw, &
                              do_pot,  &
                              eFrame, A, f, t, pot, vpair, fpair)
#endif

                         vij = vij + vpair
                         fij(1:3) = fij(1:3) + fpair(1:3)
                      endif
                   enddo inner
                enddo

                if (loop .eq. PAIR_LOOP) then
                   if (in_switching_region) then
                      swderiv = vij*dswdr/rgrp
                      fij(1) = fij(1) + swderiv*d_grp(1)
                      fij(2) = fij(2) + swderiv*d_grp(2)
                      fij(3) = fij(3) + swderiv*d_grp(3)

                      do ia=groupStartRow(i), groupStartRow(i+1)-1
                         atom1=groupListRow(ia)
                         mf = mfactRow(atom1)
#ifdef IS_MPI
                         f_Row(1,atom1) = f_Row(1,atom1) + swderiv*d_grp(1)*mf
                         f_Row(2,atom1) = f_Row(2,atom1) + swderiv*d_grp(2)*mf
                         f_Row(3,atom1) = f_Row(3,atom1) + swderiv*d_grp(3)*mf
#else
                         f(1,atom1) = f(1,atom1) + swderiv*d_grp(1)*mf
                         f(2,atom1) = f(2,atom1) + swderiv*d_grp(2)*mf
                         f(3,atom1) = f(3,atom1) + swderiv*d_grp(3)*mf
#endif
                      enddo

                      do jb=groupStartCol(j), groupStartCol(j+1)-1
                         atom2=groupListCol(jb)
                         mf = mfactCol(atom2)
#ifdef IS_MPI
                         f_Col(1,atom2) = f_Col(1,atom2) - swderiv*d_grp(1)*mf
                         f_Col(2,atom2) = f_Col(2,atom2) - swderiv*d_grp(2)*mf
                         f_Col(3,atom2) = f_Col(3,atom2) - swderiv*d_grp(3)*mf
#else
                         f(1,atom2) = f(1,atom2) - swderiv*d_grp(1)*mf
                         f(2,atom2) = f(2,atom2) - swderiv*d_grp(2)*mf
                         f(3,atom2) = f(3,atom2) - swderiv*d_grp(3)*mf
#endif
                      enddo
                   endif

                   if (do_stress) call add_stress_tensor(d_grp, fij)
                endif
             end if
          enddo

       enddo outer

       if (update_nlist) then
#ifdef IS_MPI
          point(nGroupsInRow + 1) = nlist + 1
#else 
          point(nGroups) = nlist + 1
#endif
          if (loop .eq. PREPAIR_LOOP) then
             ! we just did the neighbor list update on the first
             ! pass, so we don't need to do it
             ! again on the second pass
             update_nlist = .false.                              
          endif
       endif

       if (loop .eq. PREPAIR_LOOP) then
          call do_preforce(nlocal, pot)
       endif

    enddo

    !! Do timing
#ifdef PROFILE
    call cpu_time(forceTimeFinal)
    forceTime = forceTime + forceTimeFinal - forceTimeInitial
#endif    

#ifdef IS_MPI
    !!distribute forces

    f_temp = 0.0_dp
    call scatter(f_Row,f_temp,plan_atom_row_3d)
    do i = 1,nlocal
       f(1:3,i) = f(1:3,i) + f_temp(1:3,i)
    end do

    f_temp = 0.0_dp
    call scatter(f_Col,f_temp,plan_atom_col_3d)
    do i = 1,nlocal
       f(1:3,i) = f(1:3,i) + f_temp(1:3,i)
    end do

    if (FF_UsesDirectionalAtoms() .and. SIM_uses_DirectionalAtoms) then
       t_temp = 0.0_dp
       call scatter(t_Row,t_temp,plan_atom_row_3d)
       do i = 1,nlocal
          t(1:3,i) = t(1:3,i) + t_temp(1:3,i)
       end do
       t_temp = 0.0_dp
       call scatter(t_Col,t_temp,plan_atom_col_3d)

       do i = 1,nlocal
          t(1:3,i) = t(1:3,i) + t_temp(1:3,i)
       end do
    endif

    if (do_pot) then
       ! scatter/gather pot_row into the members of my column
       do i = 1,LR_POT_TYPES
          call scatter(pot_Row(i,:), pot_Temp(i,:), plan_atom_row)
       end do
       ! scatter/gather pot_local into all other procs
       ! add resultant to get total pot
       do i = 1, nlocal
          pot_local(1:LR_POT_TYPES) = pot_local(1:LR_POT_TYPES) &
               + pot_Temp(1:LR_POT_TYPES,i)
       enddo

       pot_Temp = 0.0_DP 
       do i = 1,LR_POT_TYPES
          call scatter(pot_Col(i,:), pot_Temp(i,:), plan_atom_col)
       end do
       do i = 1, nlocal
          pot_local(1:LR_POT_TYPES) = pot_local(1:LR_POT_TYPES)&
               + pot_Temp(1:LR_POT_TYPES,i)
       enddo

    endif
#endif

    if (FF_RequiresPostpairCalc() .and. SIM_requires_postpair_calc) then

       if (electrostaticSummationMethod == REACTION_FIELD) then

#ifdef IS_MPI
          call scatter(rf_Row,rf,plan_atom_row_3d)
          call scatter(rf_Col,rf_Temp,plan_atom_col_3d)
          do i = 1,nlocal
             rf(1:3,i) = rf(1:3,i) + rf_Temp(1:3,i)
          end do
#endif

          do i = 1, nLocal

             rfpot = 0.0_DP
#ifdef IS_MPI
             me_i = atid_row(i)
#else
             me_i = atid(i)
#endif
             iHash = InteractionHash(me_i,me_j)
             
             if ( iand(iHash, ELECTROSTATIC_PAIR).ne.0 ) then

                mu_i = getDipoleMoment(me_i)

                !! The reaction field needs to include a self contribution 
                !! to the field:
                call accumulate_self_rf(i, mu_i, eFrame)
                !! Get the reaction field contribution to the 
                !! potential and torques:
                call reaction_field_final(i, mu_i, eFrame, rfpot, t, do_pot)
#ifdef IS_MPI
                pot_local(ELECTROSTATIC_POT) = pot_local(ELECTROSTATIC_POT) + rfpot
#else
                pot(ELECTROSTATIC_POT) = pot(ELECTROSTATIC_POT) + rfpot

#endif
             endif
          enddo
       endif
    endif


#ifdef IS_MPI

    if (do_pot) then
       pot(1:LR_POT_TYPES) = pot(1:LR_POT_TYPES) &
            + pot_local(1:LR_POT_TYPES)
       !! we assume the c code will do the allreduce to get the total potential
       !! we could do it right here if we needed to...
    endif

    if (do_stress) then
       call mpi_allreduce(tau_Temp, tau, 9,mpi_double_precision,mpi_sum, &
            mpi_comm_world,mpi_err)
       call mpi_allreduce(virial_Temp, virial,1,mpi_double_precision,mpi_sum, &
            mpi_comm_world,mpi_err)
    endif

#else

    if (do_stress) then
       tau = tau_Temp
       virial = virial_Temp
    endif

#endif

  end subroutine do_force_loop

  subroutine do_pair(i, j, rijsq, d, sw, do_pot, &
       eFrame, A, f, t, pot, vpair, fpair)

    real( kind = dp ) :: vpair, sw
    real( kind = dp ), dimension(LR_POT_TYPES) :: pot
    real( kind = dp ), dimension(3) :: fpair
    real( kind = dp ), dimension(nLocal)   :: mfact
    real( kind = dp ), dimension(9,nLocal) :: eFrame
    real( kind = dp ), dimension(9,nLocal) :: A
    real( kind = dp ), dimension(3,nLocal) :: f
    real( kind = dp ), dimension(3,nLocal) :: t

    logical, intent(inout) :: do_pot
    integer, intent(in) :: i, j
    real ( kind = dp ), intent(inout) :: rijsq
    real ( kind = dp )                :: r
    real ( kind = dp ), intent(inout) :: d(3)
    integer :: me_i, me_j

    integer :: iHash

    r = sqrt(rijsq)
    vpair = 0.0d0
    fpair(1:3) = 0.0d0

#ifdef IS_MPI
    me_i = atid_row(i)
    me_j = atid_col(j)
#else
    me_i = atid(i)
    me_j = atid(j)
#endif

    iHash = InteractionHash(me_i, me_j)

    if ( iand(iHash, LJ_PAIR).ne.0 ) then
       call do_lj_pair(i, j, d, r, rijsq, sw, vpair, fpair, &
            pot(VDW_POT), f, do_pot)
    endif

    if ( iand(iHash, ELECTROSTATIC_PAIR).ne.0 ) then
       call doElectrostaticPair(i, j, d, r, rijsq, sw, vpair, fpair, &
            pot(ELECTROSTATIC_POT), eFrame, f, t, do_pot)

       if (electrostaticSummationMethod == REACTION_FIELD) then

          ! CHECK ME (RF needs to know about all electrostatic types)
          call accumulate_rf(i, j, r, eFrame, sw)
          call rf_correct_forces(i, j, d, r, eFrame, sw, f, fpair)
       endif

    endif

    if ( iand(iHash, STICKY_PAIR).ne.0 ) then
       call do_sticky_pair(i, j, d, r, rijsq, sw, vpair, fpair, &
            pot(HB_POT), A, f, t, do_pot)
    endif

    if ( iand(iHash, STICKYPOWER_PAIR).ne.0 ) then
       call do_sticky_power_pair(i, j, d, r, rijsq, sw, vpair, fpair, &
            pot(HB_POT), A, f, t, do_pot)
    endif

    if ( iand(iHash, GAYBERNE_PAIR).ne.0 ) then
       call do_gb_pair(i, j, d, r, rijsq, sw, vpair, fpair, &
            pot(VDW_POT), A, f, t, do_pot)
    endif
    
    if ( iand(iHash, GAYBERNE_LJ).ne.0 ) then
       call do_gb_lj_pair(i, j, d, r, rijsq, sw, vpair, fpair, &
            pot(VDW_POT), A, f, t, do_pot)
    endif

    if ( iand(iHash, EAM_PAIR).ne.0 ) then       
       call do_eam_pair(i, j, d, r, rijsq, sw, vpair, fpair, &
            pot(METALLIC_POT), f, do_pot)
    endif

    if ( iand(iHash, SHAPE_PAIR).ne.0 ) then       
       call do_shape_pair(i, j, d, r, rijsq, sw, vpair, fpair, &
            pot(VDW_POT), A, f, t, do_pot)
    endif

    if ( iand(iHash, SHAPE_LJ).ne.0 ) then       
       call do_shape_pair(i, j, d, r, rijsq, sw, vpair, fpair, &
            pot(VDW_POT), A, f, t, do_pot)
    endif
    
  end subroutine do_pair

  subroutine do_prepair(i, j, rijsq, d, sw, rcijsq, dc, &
       do_pot, do_stress, eFrame, A, f, t, pot)

    real( kind = dp ) :: sw
    real( kind = dp ), dimension(LR_POT_TYPES) :: pot
    real( kind = dp ), dimension(9,nLocal) :: eFrame
    real (kind=dp), dimension(9,nLocal) :: A
    real (kind=dp), dimension(3,nLocal) :: f
    real (kind=dp), dimension(3,nLocal) :: t

    logical, intent(inout) :: do_pot, do_stress
    integer, intent(in) :: i, j
    real ( kind = dp ), intent(inout)    :: rijsq, rcijsq
    real ( kind = dp )                :: r, rc
    real ( kind = dp ), intent(inout) :: d(3), dc(3)

    integer :: me_i, me_j, iHash

    r = sqrt(rijsq)

#ifdef IS_MPI   
    me_i = atid_row(i)
    me_j = atid_col(j)   
#else   
    me_i = atid(i)
    me_j = atid(j)   
#endif

    iHash = InteractionHash(me_i, me_j)

    if ( iand(iHash, EAM_PAIR).ne.0 ) then       
            call calc_EAM_prepair_rho(i, j, d, r, rijsq )
    endif
    
  end subroutine do_prepair


  subroutine do_preforce(nlocal,pot)
    integer :: nlocal
    real( kind = dp ),dimension(LR_POT_TYPES) :: pot

    if (FF_uses_EAM .and. SIM_uses_EAM) then
       call calc_EAM_preforce_Frho(nlocal,pot(METALLIC_POT))
    endif


  end subroutine do_preforce


  subroutine get_interatomic_vector(q_i, q_j, d, r_sq)

    real (kind = dp), dimension(3) :: q_i
    real (kind = dp), dimension(3) :: q_j
    real ( kind = dp ), intent(out) :: r_sq
    real( kind = dp ) :: d(3), scaled(3)
    integer i

    d(1:3) = q_j(1:3) - q_i(1:3)

    ! Wrap back into periodic box if necessary
    if ( SIM_uses_PBC ) then

       if( .not.boxIsOrthorhombic ) then
          ! calc the scaled coordinates.

          scaled = matmul(HmatInv, d)

          ! wrap the scaled coordinates

          scaled = scaled  - anint(scaled)


          ! calc the wrapped real coordinates from the wrapped scaled 
          ! coordinates

          d = matmul(Hmat,scaled)

       else
          ! calc the scaled coordinates.

          do i = 1, 3
             scaled(i) = d(i) * HmatInv(i,i)

             ! wrap the scaled coordinates

             scaled(i) = scaled(i) - anint(scaled(i))

             ! calc the wrapped real coordinates from the wrapped scaled 
             ! coordinates

             d(i) = scaled(i)*Hmat(i,i)
          enddo
       endif

    endif

    r_sq = dot_product(d,d)

  end subroutine get_interatomic_vector

  subroutine zero_work_arrays()

#ifdef IS_MPI

    q_Row = 0.0_dp
    q_Col = 0.0_dp

    q_group_Row = 0.0_dp
    q_group_Col = 0.0_dp   

    eFrame_Row = 0.0_dp
    eFrame_Col = 0.0_dp

    A_Row = 0.0_dp
    A_Col = 0.0_dp

    f_Row = 0.0_dp
    f_Col = 0.0_dp
    f_Temp = 0.0_dp

    t_Row = 0.0_dp
    t_Col = 0.0_dp
    t_Temp = 0.0_dp

    pot_Row = 0.0_dp
    pot_Col = 0.0_dp
    pot_Temp = 0.0_dp

    rf_Row = 0.0_dp
    rf_Col = 0.0_dp
    rf_Temp = 0.0_dp

#endif

    if (FF_uses_EAM .and. SIM_uses_EAM) then
       call clean_EAM()
    endif

    rf = 0.0_dp
    tau_Temp = 0.0_dp
    virial_Temp = 0.0_dp
  end subroutine zero_work_arrays

  function skipThisPair(atom1, atom2) result(skip_it)
    integer, intent(in) :: atom1
    integer, intent(in), optional :: atom2
    logical :: skip_it
    integer :: unique_id_1, unique_id_2
    integer :: me_i,me_j
    integer :: i

    skip_it = .false. 

    !! there are a number of reasons to skip a pair or a particle
    !! mostly we do this to exclude atoms who are involved in short
    !! range interactions (bonds, bends, torsions), but we also need 
    !! to exclude some overcounted interactions that result from
    !! the parallel decomposition

#ifdef IS_MPI
    !! in MPI, we have to look up the unique IDs for each atom
    unique_id_1 = AtomRowToGlobal(atom1)
#else
    !! in the normal loop, the atom numbers are unique
    unique_id_1 = atom1
#endif

    !! We were called with only one atom, so just check the global exclude
    !! list for this atom
    if (.not. present(atom2)) then
       do i = 1, nExcludes_global
          if (excludesGlobal(i) == unique_id_1) then
             skip_it = .true.
             return
          end if
       end do
       return 
    end if

#ifdef IS_MPI
    unique_id_2 = AtomColToGlobal(atom2)
#else
    unique_id_2 = atom2
#endif

#ifdef IS_MPI
    !! this situation should only arise in MPI simulations
    if (unique_id_1 == unique_id_2) then
       skip_it = .true.
       return
    end if

    !! this prevents us from doing the pair on multiple processors
    if (unique_id_1 < unique_id_2) then
       if (mod(unique_id_1 + unique_id_2,2) == 0) then
          skip_it = .true.
          return
       endif
    else                
       if (mod(unique_id_1 + unique_id_2,2) == 1) then 
          skip_it = .true.
          return
       endif
    endif
#endif

    !! the rest of these situations can happen in all simulations:
    do i = 1, nExcludes_global       
       if ((excludesGlobal(i) == unique_id_1) .or. &
            (excludesGlobal(i) == unique_id_2)) then
          skip_it = .true.
          return
       endif
    enddo

    do i = 1, nSkipsForAtom(atom1)
       if (skipsForAtom(atom1, i) .eq. unique_id_2) then
          skip_it = .true.
          return
       endif
    end do

    return
  end function skipThisPair

  function FF_UsesDirectionalAtoms() result(doesit)
    logical :: doesit
    doesit = FF_uses_DirectionalAtoms
  end function FF_UsesDirectionalAtoms

  function FF_RequiresPrepairCalc() result(doesit)
    logical :: doesit
    doesit = FF_uses_EAM
  end function FF_RequiresPrepairCalc

  function FF_RequiresPostpairCalc() result(doesit)
    logical :: doesit
    if (electrostaticSummationMethod == REACTION_FIELD) doesit = .true.
  end function FF_RequiresPostpairCalc

#ifdef PROFILE
  function getforcetime() result(totalforcetime)
    real(kind=dp) :: totalforcetime
    totalforcetime = forcetime
  end function getforcetime
#endif

  !! This cleans componets of force arrays belonging only to fortran

  subroutine add_stress_tensor(dpair, fpair)

    real( kind = dp ), dimension(3), intent(in) :: dpair, fpair

    ! because the d vector is the rj - ri vector, and
    ! because fx, fy, fz are the force on atom i, we need a
    ! negative sign here:  

    tau_Temp(1) = tau_Temp(1) - dpair(1) * fpair(1)
    tau_Temp(2) = tau_Temp(2) - dpair(1) * fpair(2)
    tau_Temp(3) = tau_Temp(3) - dpair(1) * fpair(3)
    tau_Temp(4) = tau_Temp(4) - dpair(2) * fpair(1)
    tau_Temp(5) = tau_Temp(5) - dpair(2) * fpair(2)
    tau_Temp(6) = tau_Temp(6) - dpair(2) * fpair(3)
    tau_Temp(7) = tau_Temp(7) - dpair(3) * fpair(1)
    tau_Temp(8) = tau_Temp(8) - dpair(3) * fpair(2)
    tau_Temp(9) = tau_Temp(9) - dpair(3) * fpair(3)

    virial_Temp = virial_Temp + &
         (tau_Temp(1) + tau_Temp(5) + tau_Temp(9))

  end subroutine add_stress_tensor

end module doForces
