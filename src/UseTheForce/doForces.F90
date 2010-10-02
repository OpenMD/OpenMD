!!
!! Copyright (c) 2005, 2009 The University of Notre Dame. All Rights Reserved.
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

!! doForces.F90
!! module doForces
!! Calculates Long Range forces.

!! @author Charles F. Vardeman II
!! @author Matthew Meineke
!! @version $Id$, $Date$, $Name: not supported by cvs2svn $, $Revision$


module doForces
  use force_globals
  use fForceOptions
  use simulation
  use definitions
  use atype_module
  use switcheroo
  use neighborLists  
  use shapes
  use vector_class
  use MetalNonMetal
  use status
  use ISO_C_BINDING

#ifdef IS_MPI
  use mpiSimulation
#endif

  implicit none
  PRIVATE

  real(kind=dp), external :: getSigma
  real(kind=dp), external :: getEpsilon
  real(kind=dp), external :: getEAMcut
  real(kind=dp), external :: getGayBerneCut
  real(kind=dp), external :: getStickyCut
  real(kind=dp), external :: getSCCut

  
#define __FORTRAN90
#include "UseTheForce/fCutoffPolicy.h"
#include "UseTheForce/DarkSide/fInteractionMap.h"
#include "UseTheForce/DarkSide/fElectrostaticSummationMethod.h"

  INTEGER, PARAMETER:: PREPAIR_LOOP = 1
  INTEGER, PARAMETER:: PAIR_LOOP    = 2

  logical, save :: haveNeighborList = .false.
  logical, save :: haveSIMvariables = .false.
  logical, save :: haveSaneForceField = .false.
  logical, save :: haveGtypeCutoffMap = .false.
  logical, save :: haveDefaultCutoffs = .false.
  logical, save :: haveSkinThickness = .false.
  logical, save :: haveElectrostaticSummationMethod = .false.
  logical, save :: haveCutoffPolicy = .false.
  logical, save :: VisitCutoffsAfterComputing = .false.

  logical, save :: FF_uses_DirectionalAtoms
  logical, save :: FF_uses_Dipoles
  logical, save :: FF_uses_GayBerne
  logical, save :: FF_uses_EAM
  logical, save :: FF_uses_SC
  logical, save :: FF_uses_MNM
 

  logical, save :: SIM_uses_DirectionalAtoms
  logical, save :: SIM_uses_EAM
  logical, save :: SIM_uses_SC
  logical, save :: SIM_uses_MNM
  logical, save :: SIM_requires_postpair_calc
  logical, save :: SIM_requires_prepair_calc
  logical, save :: SIM_uses_PBC
  logical, save :: SIM_uses_AtomicVirial

  integer, save :: electrostaticSummationMethod
  integer, save :: cutoffPolicy = TRADITIONAL_CUTOFF_POLICY

  real(kind=dp), save :: defaultRcut, defaultRsw, largestRcut
  real(kind=dp), save :: skinThickness
  logical, save :: defaultDoShiftPot
  logical, save :: defaultDoShiftFrc

  public :: init_FF
  public :: setCutoffs
  public :: cWasLame
  public :: setElectrostaticMethod
  public :: setCutoffPolicy
  public :: setSkinThickness
  public :: do_force_loop

#ifdef PROFILE
  public :: getforcetime
  real, save :: forceTime = 0
  real :: forceTimeInitial, forceTimeFinal
  integer :: nLoops
#endif
  
  !! Variables for cutoff mapping and interaction mapping
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
   
contains

  subroutine createGtypeCutoffMap()

    logical :: i_is_LJ
    logical :: i_is_Elect
    logical :: i_is_Sticky
    logical :: i_is_StickyP
    logical :: i_is_GB
    logical :: i_is_EAM
    logical :: i_is_Shape
    logical :: i_is_SC
    logical :: GtypeFound

    integer :: myStatus, nAtypes,  i, j, istart, iend, jstart, jend
    integer :: n_in_i, me_i, ia, g, atom1, ja, n_in_j,me_j
    integer :: nGroupsInRow
    integer :: nGroupsInCol
    integer :: nGroupTypesRow,nGroupTypesCol
    real(kind=dp):: thisSigma, bigSigma, thisRcut, tradRcut, tol
    real(kind=dp) :: biggestAtypeCutoff
    integer :: c_ident_i

#ifdef IS_MPI
    nGroupsInRow = getNgroupsInRow(plan_group_row)
    nGroupsInCol = getNgroupsInCol(plan_group_col)
#endif
    nAtypes = getSize(atypes)
! Set all of the initial cutoffs to zero.
    atypeMaxCutoff = 0.0_dp
    biggestAtypeCutoff = 0.0_dp
    do i = 1, nAtypes
       if (SimHasAtype(i)) then     
          call getElementProperty(atypes, i, "is_LennardJones", i_is_LJ)
          call getElementProperty(atypes, i, "is_Electrostatic", i_is_Elect)
          call getElementProperty(atypes, i, "is_Sticky", i_is_Sticky)
          call getElementProperty(atypes, i, "is_StickyPower", i_is_StickyP)
          call getElementProperty(atypes, i, "is_GayBerne", i_is_GB)
          call getElementProperty(atypes, i, "is_EAM", i_is_EAM)
          call getElementProperty(atypes, i, "is_Shape", i_is_Shape)
          call getElementProperty(atypes, i, "is_SC", i_is_SC)
          call getElementProperty(atypes, i, "c_ident", c_ident_i)
 
          if (haveDefaultCutoffs) then
             atypeMaxCutoff(i) = defaultRcut
          else
             if (i_is_LJ) then           
                thisRcut = getSigma(c_ident_i)  * 2.5_dp
                if (thisRCut .gt. atypeMaxCutoff(i)) atypeMaxCutoff(i) = thisRCut
             endif
             if (i_is_Elect) then
                thisRcut = defaultRcut
                if (thisRCut .gt. atypeMaxCutoff(i)) atypeMaxCutoff(i) = thisRCut
             endif
             if (i_is_Sticky) then
                thisRcut = getStickyCut(c_ident_i)
                if (thisRCut .gt. atypeMaxCutoff(i)) atypeMaxCutoff(i) = thisRCut
             endif
             if (i_is_StickyP) then
                thisRcut = getStickyCut(c_ident_i)
                if (thisRCut .gt. atypeMaxCutoff(i)) atypeMaxCutoff(i) = thisRCut
             endif
             if (i_is_GB) then
                thisRcut = getGayBerneCut(c_ident_i)
                if (thisRCut .gt. atypeMaxCutoff(i)) atypeMaxCutoff(i) = thisRCut
             endif
             if (i_is_EAM) then
                thisRcut = getEAMCut(c_ident_i)
                if (thisRCut .gt. atypeMaxCutoff(i)) atypeMaxCutoff(i) = thisRCut
             endif
             if (i_is_Shape) then
                thisRcut = getShapeCut(i)
                if (thisRCut .gt. atypeMaxCutoff(i)) atypeMaxCutoff(i) = thisRCut
             endif
             if (i_is_SC) then
                thisRcut = getSCCut(c_ident_i)
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

    if(.not.associated(groupMaxCutoffCol)) then
       allocate(groupMaxCutoffCol(jend))
    else
       deallocate(groupMaxCutoffCol)
       allocate(groupMaxCutoffCol(jend))
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
    
    tol = 1.0e-6_dp
    nGroupTypesRow = 0
    nGroupTypesCol = 0
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
          
          if (thisRcut.gt.largestRcut) largestRcut = thisRcut

          gtypeCutoffMap(i,j)%rcutsq = thisRcut*thisRcut

          if (.not.haveSkinThickness) then
             skinThickness = 1.0_dp
          endif

          gtypeCutoffMap(i,j)%rlistsq = (thisRcut + skinThickness)**2

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

   subroutine setCutoffs(defRcut, defRsw, defSP, defSF)

     real(kind=dp),intent(in) :: defRcut, defRsw
     integer, intent(in) :: defSP, defSF
     character(len = statusMsgSize) :: errMsg
     integer :: localError

     defaultRcut = defRcut
     defaultRsw = defRsw
    
     if (defSP .ne. 0) then  
        defaultDoShiftPot = .true.
     else
        defaultDoShiftPot = .false.
     endif
     if (defSF .ne. 0) then  
        defaultDoShiftFrc = .true.
     else
        defaultDoShiftFrc = .false.
     endif

     if (abs(defaultRcut-defaultRsw) .lt. 0.0001) then
        if (defaultDoShiftFrc) then
           write(errMsg, *) &
                'cutoffRadius and switchingRadius are set to the', newline &
                // tab, 'same value.  OpenMD will use shifted force', newline &
                // tab, 'potentials instead of switching functions.'
           
           call handleInfo("setCutoffs", errMsg)
        else
           write(errMsg, *) &
                'cutoffRadius and switchingRadius are set to the', newline &
                // tab, 'same value.  OpenMD will use shifted', newline &
                // tab, 'potentials instead of switching functions.'
           
           call handleInfo("setCutoffs", errMsg)
           
           defaultDoShiftPot = .true.
        endif
                
     endif
     
     localError = 0
     call setLJDefaultCutoff( defaultRcut, defaultDoShiftPot, &
          defaultDoShiftFrc )
     call setElectrostaticCutoffRadius( defaultRcut, defaultRsw )
     call setCutoffEAM( defaultRcut )
     call setCutoffSC( defaultRcut )
     call setMnMDefaultCutoff( defaultRcut, defaultDoShiftPot, &
          defaultDoShiftFrc )
     call set_switch(defaultRsw, defaultRcut)
     call setHmatDangerousRcutValue(defaultRcut)
         
     haveDefaultCutoffs = .true.
     haveGtypeCutoffMap = .false.

   end subroutine setCutoffs

   subroutine cWasLame()
     
     VisitCutoffsAfterComputing = .true.
     return
     
   end subroutine cWasLame
   
   subroutine setCutoffPolicy(cutPolicy)
     
     integer, intent(in) :: cutPolicy
     
     cutoffPolicy = cutPolicy
     haveCutoffPolicy = .true.
     haveGtypeCutoffMap = .false.
     
   end subroutine setCutoffPolicy
     

   subroutine setElectrostaticMethod( thisESM )

     integer, intent(in) :: thisESM

     electrostaticSummationMethod = thisESM
     haveElectrostaticSummationMethod = .true.
    
   end subroutine setElectrostaticMethod

   subroutine setSkinThickness( thisSkin )
     
     real(kind=dp), intent(in) :: thisSkin
     
     skinThickness = thisSkin
     haveSkinThickness = .true.    
     haveGtypeCutoffMap = .false.
     
   end subroutine setSkinThickness
      
   subroutine setSimVariables()
     SIM_uses_DirectionalAtoms = SimUsesDirectionalAtoms()
     SIM_uses_EAM = SimUsesEAM()
     SIM_requires_postpair_calc = SimRequiresPostpairCalc()
     SIM_requires_prepair_calc = SimRequiresPrepairCalc()
     SIM_uses_PBC = SimUsesPBC()
     SIM_uses_SC = SimUsesSC()
     SIM_uses_AtomicVirial = SimUsesAtomicVirial()

     haveSIMvariables = .true.
     
     return
   end subroutine setSimVariables

  subroutine doReadyCheck(error)
    integer, intent(out) :: error
    integer :: myStatus

    error = 0

    if (.not. haveGtypeCutoffMap) then        
       call createGtypeCutoffMap()       
    endif

    if (VisitCutoffsAfterComputing) then
       call set_switch(largestRcut, largestRcut)       
       call setHmatDangerousRcutValue(largestRcut)
       call setCutoffEAM(largestRcut)
       call setCutoffSC(largestRcut)
       VisitCutoffsAfterComputing = .false.
    endif

    if (.not. haveSIMvariables) then
       call setSimVariables()
    endif

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


  subroutine init_FF(thisStat)

    integer, intent(out) :: thisStat   
    integer :: my_status, nMatches
    integer, pointer :: MatchList(:) => null()

    !! assume things are copacetic, unless they aren't
    thisStat = 0

    !! init_FF is called *after* all of the atom types have been 
    !! defined in atype_module using the new_atype subroutine.
    !!
    !! this will scan through the known atypes and figure out what
    !! interactions are used by the force field.     

    FF_uses_DirectionalAtoms = .false.
    FF_uses_Dipoles = .false.
    FF_uses_GayBerne = .false.
    FF_uses_EAM = .false.
    FF_uses_SC = .false.

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

    call getMatchingElementList(atypes, "is_SC", .true., nMatches, MatchList)
    if (nMatches .gt. 0) FF_uses_SC = .true.


    haveSaneForceField = .true.


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
  subroutine do_force_loop(q, q_group, A, eFrame, f, t, tau, pot, particle_pot, &
       error)
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
    real( kind = dp ), dimension(nLocal) :: particle_pot

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
    real( kind = DP ) :: ratmsq, rgrpsq, rgrp, rag, vpair, vij
    real( kind = DP ) :: sw, dswdr, swderiv, mf
    real( kind = DP ) :: rVal
    real(kind=dp),dimension(3) :: d_atm, d_grp, fpair, fij, fg, dag
    real(kind=dp) :: rfpot, mu_i
    real(kind=dp):: rCut
    integer :: me_i, me_j, n_in_i, n_in_j, iG, j1
    logical :: is_dp_i
    integer :: neighborListSize
    integer :: listerror, error
    integer :: localError
    integer :: propPack_i, propPack_j
    integer :: loopStart, loopEnd, loop
    integer :: i1, topoDist

    real(kind=dp) :: skch

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
          call checkNeighborList(nGroupsInRow, q_group_row, skinThickness, &
               update_nlist)
#else
          call checkNeighborList(nGroups, q_group, skinThickness, &
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
                 
                if (rgrpsq < gtypeCutoffMap(groupToGtypeRow(i),groupToGtypeCol(j))%rCutsq) then

                   rCut = gtypeCutoffMap(groupToGtypeRow(i),groupToGtypeCol(j))%rCut
                   if (loop .eq. PAIR_LOOP) then
                      vij = 0.0_dp
                      fij(1) = 0.0_dp
                      fij(2) = 0.0_dp
                      fij(3) = 0.0_dp
                   endif
                   
                   call get_switch(rgrpsq, sw, dswdr,rgrp, in_switching_region)
                   
                   n_in_j = groupStartCol(j+1) - groupStartCol(j)
                   
                   do ia = groupStartRow(i), groupStartRow(i+1)-1
                      
                      atom1 = groupListRow(ia)
                      
                      inner: do jb = groupStartCol(j), groupStartCol(j+1)-1
                         
                         atom2 = groupListCol(jb)
                         
                         if (skipThisPair(atom1, atom2))  cycle inner
                         
                         if ((n_in_i .eq. 1).and.(n_in_j .eq. 1)) then
                            d_atm(1) = d_grp(1)
                            d_atm(2) = d_grp(2)
                            d_atm(3) = d_grp(3)
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

                         topoDist = getTopoDistance(atom1, atom2)

                         if (loop .eq. PREPAIR_LOOP) then
#ifdef IS_MPI                      
                            call do_prepair(atom1, atom2, ratmsq, d_atm, sw, &
                                 rgrpsq, d_grp, rCut, &
                                 eFrame, A, f, t, pot_local)
#else
                            call do_prepair(atom1, atom2, ratmsq, d_atm, sw, &
                                 rgrpsq, d_grp, rCut, &
                                 eFrame, A, f, t, pot)
#endif                                               
                         else
#ifdef IS_MPI                      
                            call do_pair(atom1, atom2, ratmsq, d_atm, sw, &
                                 eFrame, A, f, t, pot_local, particle_pot, vpair, &
                                 fpair, d_grp, rgrp, rCut, topoDist)
                            ! particle_pot will be accumulated from row & column
                            ! arrays later
#else
                            call do_pair(atom1, atom2, ratmsq, d_atm, sw, &
                                 eFrame, A, f, t, pot, particle_pot, vpair, &
                                 fpair, d_grp, rgrp, rCut, topoDist)
#endif
                            vij = vij + vpair
                            fij(1) = fij(1) + fpair(1)
                            fij(2) = fij(2) + fpair(2)
                            fij(3) = fij(3) + fpair(3)
                            call add_stress_tensor(d_atm, fpair, tau)
                         endif
                      enddo inner
                   enddo

                   if (loop .eq. PAIR_LOOP) then
                      if (in_switching_region) then
                         swderiv = vij*dswdr/rgrp
                         fg = swderiv*d_grp
 
                         fij(1) = fij(1) + fg(1)
                         fij(2) = fij(2) + fg(2)
                         fij(3) = fij(3) + fg(3)
                         
                         if ((n_in_i .eq. 1).and.(n_in_j .eq. 1)) then
                            call add_stress_tensor(d_atm, fg, tau)
                         endif
                         
                         do ia=groupStartRow(i), groupStartRow(i+1)-1
                            atom1=groupListRow(ia)
                            mf = mfactRow(atom1)
                            ! fg is the force on atom ia due to cutoff group's
                            ! presence in switching region
                            fg = swderiv*d_grp*mf
#ifdef IS_MPI
                            f_Row(1,atom1) = f_Row(1,atom1) + fg(1)
                            f_Row(2,atom1) = f_Row(2,atom1) + fg(2)
                            f_Row(3,atom1) = f_Row(3,atom1) + fg(3)
#else
                            f(1,atom1) = f(1,atom1) + fg(1)
                            f(2,atom1) = f(2,atom1) + fg(2)
                            f(3,atom1) = f(3,atom1) + fg(3)
#endif
                            if (n_in_i .gt. 1) then
                               if (SIM_uses_AtomicVirial) then
                                  ! find the distance between the atom
                                  ! and the center of the cutoff group:
#ifdef IS_MPI
                                  call get_interatomic_vector(q_Row(:,atom1), &
                                       q_group_Row(:,i), dag, rag)
#else
                                  call get_interatomic_vector(q(:,atom1), &
                                       q_group(:,i), dag, rag)
#endif
                                  call add_stress_tensor(dag,fg,tau)
                               endif
                            endif
                         enddo
                         
                         do jb=groupStartCol(j), groupStartCol(j+1)-1
                            atom2=groupListCol(jb)
                            mf = mfactCol(atom2)
                            ! fg is the force on atom jb due to cutoff group's
                            ! presence in switching region
                            fg = -swderiv*d_grp*mf
#ifdef IS_MPI
                            f_Col(1,atom2) = f_Col(1,atom2) + fg(1)
                            f_Col(2,atom2) = f_Col(2,atom2) + fg(2)
                            f_Col(3,atom2) = f_Col(3,atom2) + fg(3)
#else
                            f(1,atom2) = f(1,atom2) + fg(1)
                            f(2,atom2) = f(2,atom2) + fg(2)
                            f(3,atom2) = f(3,atom2) + fg(3)
#endif
                            if (n_in_j .gt. 1) then
                               if (SIM_uses_AtomicVirial) then
                                  ! find the distance between the atom
                                  ! and the center of the cutoff group:
#ifdef IS_MPI
                                  call get_interatomic_vector(q_Col(:,atom2), &
                                       q_group_Col(:,j), dag, rag)
#else
                                  call get_interatomic_vector(q(:,atom2), &
                                       q_group(:,j), dag, rag)
#endif
                                  call add_stress_tensor(dag,fg,tau)
                               endif
                            endif
                         enddo
                      endif
                      !if (.not.SIM_uses_AtomicVirial) then
                      !   call add_stress_tensor(d_grp, fij, tau) 
                      !endif
                   endif
                endif
             endif
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
#ifdef IS_MPI
          call do_preforce(nlocal, pot_local, particle_pot)
#else
          call do_preforce(nlocal, pot, particle_pot)
#endif
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
    
    do i = 1,LR_POT_TYPES
       particle_pot(1:nlocal) = particle_pot(1:nlocal) + pot_Temp(i,1:nlocal)
    enddo
    
    pot_Temp = 0.0_DP 
    
    do i = 1,LR_POT_TYPES
       call scatter(pot_Col(i,:), pot_Temp(i,:), plan_atom_col)
    end do
    
    do i = 1, nlocal
       pot_local(1:LR_POT_TYPES) = pot_local(1:LR_POT_TYPES)&
            + pot_Temp(1:LR_POT_TYPES,i)
    enddo
    
    do i = 1,LR_POT_TYPES
       particle_pot(1:nlocal) = particle_pot(1:nlocal) + pot_Temp(i,1:nlocal)
    enddo
    
    ppot_Temp = 0.0_DP
    
    call scatter(ppot_Row(:), ppot_Temp(:), plan_atom_row)
    do i = 1, nlocal
       particle_pot(i) = particle_pot(i) + ppot_Temp(i)
    enddo
    
    ppot_Temp = 0.0_DP
    
    call scatter(ppot_Col(:), ppot_Temp(:), plan_atom_col)
    do i = 1, nlocal
       particle_pot(i) = particle_pot(i) + ppot_Temp(i)
    enddo
   
#endif

    if (SIM_requires_postpair_calc) then
       do i = 1, nlocal             
          
          ! we loop only over the local atoms, so we don't need row and column
          ! lookups for the types

          me_i = atid(i)
          
          ! is the atom electrostatic?  See if it would have an 
          ! electrostatic interaction with itself
          iHash = InteractionHash(me_i,me_i)

          if ( iand(iHash, ELECTROSTATIC_PAIR).ne.0 ) then

             ! loop over the excludes to accumulate charge in the
             ! cutoff sphere that we've left out of the normal pair loop
             skch = 0.0_dp
                          
             do i1 = 1, nSkipsForLocalAtom(i)                
                j = skipsForLocalAtom(i, i1)                
                me_j = atid(j)
                jHash = InteractionHash(me_i,me_j)
                if ( iand(jHash, ELECTROSTATIC_PAIR).ne.0 ) then
                   skch = skch + getCharge(me_j)
                endif
             enddo

#ifdef IS_MPI
             call self_self(i, eFrame, skch, pot_local(ELECTROSTATIC_POT), t)
#else
             call self_self(i, eFrame, skch, pot(ELECTROSTATIC_POT), t)
#endif
          endif
  
          
          if (electrostaticSummationMethod.eq.REACTION_FIELD) then
             
             ! loop over the excludes to accumulate RF stuff we've
             ! left out of the normal pair loop
             
             do i1 = 1, nSkipsForLocalAtom(i)
                j = skipsForLocalAtom(i, i1)
                
                ! prevent overcounting of the skips
                if (i.lt.j) then
                   call get_interatomic_vector(q(:,i), q(:,j), d_atm, ratmsq)
                   rVal = sqrt(ratmsq)
                   call get_switch(ratmsq, sw, dswdr, rVal,in_switching_region)
#ifdef IS_MPI
                   call rf_self_excludes(i, j, sw, 1.0_dp, eFrame, d_atm, rVal, &
                        vpair, pot_local(ELECTROSTATIC_POT), f, t)
#else
                   call rf_self_excludes(i, j, sw, 1.0_dp, eFrame, d_atm, rVal, &
                        vpair, pot(ELECTROSTATIC_POT), f, t)
#endif
                endif
             enddo
          endif
       enddo
    endif

#ifdef IS_MPI
#ifdef SINGLE_PRECISION
    call mpi_allreduce(pot_local, pot, LR_POT_TYPES,mpi_real,mpi_sum, &
         mpi_comm_world,mpi_err)            
#else
    call mpi_allreduce(pot_local, pot, LR_POT_TYPES,mpi_double_precision, &
         mpi_sum, mpi_comm_world,mpi_err)            
#endif
#endif

  end subroutine do_force_loop

  subroutine do_pair(i, j, rijsq, d, sw, &
       eFrame, A, f, t, pot, particle_pot, vpair, &
       fpair, d_grp, r_grp, rCut, topoDist)

    real( kind = dp ) :: vpair, sw
    real( kind = dp ), dimension(LR_POT_TYPES) :: pot, pairpot
    real( kind = dp ), dimension(nLocal) :: particle_pot
    real( kind = dp ), dimension(3) :: fpair
    real( kind = dp ), dimension(nLocal)   :: mfact
    real( kind = dp ), dimension(9,nLocal) :: eFrame
    real( kind = dp ), dimension(9,nLocal) :: A
    real( kind = dp ), dimension(3,nLocal) :: f
    real( kind = dp ), dimension(3,nLocal) :: t

    integer, intent(in) :: i, j
    real ( kind = dp ), intent(inout) :: rijsq
    real ( kind = dp ), intent(inout) :: r_grp
    real ( kind = dp ), intent(inout) :: d(3)
    real ( kind = dp ), intent(inout) :: d_grp(3)
    real ( kind = dp ), intent(inout) :: rCut 
    integer, intent(inout) :: topoDist
    real ( kind = dp ) :: r, pair_pot, vdwMult, electroMult
    real ( kind = dp ) :: a_k, b_k, c_k, d_k, dx

    real( kind = dp), dimension(3) :: f1, t1, t2
    real( kind = dp), dimension(9) :: A1, A2, eF1, eF2
    real( kind = dp) :: dfrhodrho_i, dfrhodrho_j
    real( kind = dp) :: rho_i, rho_j
    real( kind = dp) :: fshift_i, fshift_j 
    integer :: id1, id2, idx
    integer :: k
    integer :: c_ident_i, c_ident_j

    integer :: iHash

    r = sqrt(rijsq)
    
    vpair = 0.0_dp
    fpair(1:3) = 0.0_dp

    p_vdw = 0.0
    p_elect = 0.0
    p_hb = 0.0
    p_met = 0.0

    f1(1:3) = 0.0
    t1(1:3) = 0.0
    t2(1:3) = 0.0

#ifdef IS_MPI
    c_ident_i = c_idents_row(i)
    c_ident_j = c_idents_col(j)

    do idx = 1, 9
       A1(idx) = A_Row(idx, i)
       A2(idx) = A_Col(idx, j)
       eF1(idx) = eFrame_Row(idx, i)
       eF2(idx) = eFrame_Col(idx, j)
    enddo

    dfrhodrho_i = dfrhodrho_row(i)
    dfrhodrho_j = dfrhodrho_col(j)
    rho_i = rho_row(i)
    rho_j = rho_col(j)
       
#else
    c_ident_i = c_idents_local(i)
    c_ident_j = c_idents_local(j)

    do idx = 1, 9
       A1(idx) = A(idx, i)
       A2(idx) = A(idx, j)
       eF1(idx) = eFrame(idx, i)
       eF2(idx) = eFrame(idx, j)       
    enddo

    dfrhodrho_i = dfrhodrho(i)
    dfrhodrho_j = dfrhodrho(j)
    rho_i = rho(i)
    rho_j = rho(j)
    
#endif
    
    vdwMult = vdwScale(topoDist)
    electroMult = electrostaticScale(topoDist)

    call doPairInteraction(c_ident_i, c_ident_j, d, r, rijsq, sw, vpair, &
         vdwMult, electroMult, A1, A2, eF1, eF2,  &
         pairpot, f1, t1, t2, &
         rho_i, rho_j, dfrhodrho_i, dfrhodrho_j, fshift_i, fshift_j)
    
#ifdef IS_MPI
    id1 = AtomRowToGlobal(i)
    id2 = AtomColToGlobal(j)

    pot_row(VDW_POT,i) = pot_row(VDW_POT,i) + 0.5*pairpot(VDW_POT)
    pot_col(VDW_POT,j) = pot_col(VDW_POT,j) + 0.5*pairpot(VDW_POT)
    pot_row(ELECTROSTATIC_POT,i) = pot_row(ELECTROSTATIC_POT,i) + 0.5*pairpot(ELECTROSTATIC_POT)
    pot_col(ELECTROSTATIC_POT,j) = pot_col(ELECTROSTATIC_POT,j) + 0.5*pairpot(ELECTROSTATIC_POT)
    pot_row(HB_POT,i) = pot_row(HB_POT,i) + 0.5*pairpot(HB_POT)
    pot_col(HB_POT,j) = pot_col(HB_POT,j) + 0.5*pairpot(HB_POT)
    pot_Row(METALLIC_POT,i) = pot_Row(METALLIC_POT,i) + 0.5*pairpot(METALLIC_POT)
    pot_Col(METALLIC_POT,j) = pot_Col(METALLIC_POT,j) + 0.5*p(METALLIC_POT)

    do idx = 1, 3
       f_Row(idx,i) = f_Row(idx,i) + f1(idx)
       f_Col(idx,j) = f_Col(idx,j) - f1(idx)
    
       t_Row(idx,i) = t_Row(idx,i) + t1(idx)
       t_Col(idx,j) = t_Col(idx,j) + t2(idx)
    enddo
       ! particle_pot is the difference between the full potential 
       ! and the full potential without the presence of a particular
       ! particle (atom1).
       !
       ! This reduces the density at other particle locations, so
       ! we need to recompute the density at atom2 assuming atom1
       ! didn't contribute.  This then requires recomputing the
       ! density functional for atom2 as well.
       !
       ! Most of the particle_pot heavy lifting comes from the
       ! pair interaction, and will be handled by vpair. Parallel version.
       
    if ( (iand(iHash, EAM_PAIR).ne.0) .or. (iand(iHash, SC_PAIR).ne.0)  ) then
       ppot_row(i) = ppot_row(i) - frho_row(j) + fshift_j
       ppot_col(j) = ppot_col(j) - frho_col(i) + fshift_i
    end if
    
#else
    id1 = i
    id2 = j

    pot(VDW_POT) = pot(VDW_POT) + pairpot(VDW_POT)
    pot(ELECTROSTATIC_POT) = pot(ELECTROSTATIC_POT) + pairpot(ELECTROSTATIC_POT)
    pot(HB_POT) = pot(HB_POT) + pairpot(HB_POT)
    pot(METALLIC_POT) = pot(METALLIC_POT) + pairpot(METALLIC_POT)

    do idx = 1, 3
       f(idx,i) = f(idx,i) + f1(idx)
       f(idx,j) = f(idx,j) - f1(idx)

       t(idx,i) = t(idx,i) + t1(idx)
       t(idx,j) = t(idx,j) + t2(idx)
    enddo
       ! particle_pot is the difference between the full potential 
       ! and the full potential without the presence of a particular
       ! particle (atom1).
       !
       ! This reduces the density at other particle locations, so
       ! we need to recompute the density at atom2 assuming atom1
       ! didn't contribute.  This then requires recomputing the
       ! density functional for atom2 as well.
       !
       ! Most of the particle_pot heavy lifting comes from the
       ! pair interaction, and will be handled by vpair. NonParallel version.

    if ( (iand(iHash, EAM_PAIR).ne.0) .or. (iand(iHash, SC_PAIR).ne.0)  ) then
       particle_pot(i) = particle_pot(i) - frho(j) + fshift_j
       particle_pot(j) = particle_pot(j) - frho(i) + fshift_i
    end if


#endif
    
    if (molMembershipList(id1) .ne. molMembershipList(id2)) then
       
       fpair(1) = fpair(1) + f1(1)
       fpair(2) = fpair(2) + f1(2)
       fpair(3) = fpair(3) + f1(3)
       
    endif
  end subroutine do_pair

  subroutine do_prepair(i, j, rijsq, d, sw, rcijsq, dc, rCut, &
       eFrame, A, f, t, pot)
    
    real( kind = dp ) :: sw
    real( kind = dp ), dimension(LR_POT_TYPES) :: pot
    real( kind = dp ), dimension(9,nLocal) :: eFrame
    real (kind=dp), dimension(9,nLocal) :: A
    real (kind=dp), dimension(3,nLocal) :: f
    real (kind=dp), dimension(3,nLocal) :: t
    
    integer, intent(in) :: i, j
    real ( kind = dp ), intent(inout)    :: rijsq, rCut
    real ( kind = dp )                :: r
    real ( kind = dp ), intent(inout) :: d(3), dc(3)
    real ( kind = dp ) :: rho_i_at_j, rho_j_at_i
    integer ::  c_ident_i, c_ident_j
    
    r = sqrt(rijsq)
    
#ifdef IS_MPI   
    c_ident_i = c_idents_row(i)
    c_ident_j = c_idents_col(j)    
#else   
    c_ident_i = c_idents_local(i)
    c_ident_j = c_idents_local(j)    
#endif
    rho_i_at_j = 0.0_dp
    rho_j_at_i = 0.0_dp
    
    call doPrepairInteraction(c_ident_i, c_ident_j, r, &
         rho_i_at_j, rho_j_at_i)
    
#ifdef IS_MPI
    rho_col(j) = rho_col(j) + rho_i_at_j
    rho_row(i) = rho_row(i) + rho_j_at_i
#else
    rho(j) = rho(j) + rho_i_at_j
    rho(i) = rho(i) + rho_j_at_i
#endif       
    
  end subroutine do_prepair


  subroutine do_preforce(nlocal, pot, particle_pot)
    integer :: nlocal
    real( kind = dp ),dimension(LR_POT_TYPES) :: pot
    real( kind = dp ),dimension(nlocal) :: particle_pot
    integer :: sc_err = 0
    integer :: atid1, atom, c_ident1
    
    if ((FF_uses_EAM .and. SIM_uses_EAM) .or. (FF_uses_SC .and. SIM_uses_SC)) then
       
#ifdef IS_MPI
       call scatter(rho_row,rho,plan_atom_row,sc_err)
       if (sc_err /= 0 ) then
          call handleError("do_preforce()", "Error scattering rho_row into rho")
       endif
       call scatter(rho_col,rho_tmp,plan_atom_col,sc_err)
       if (sc_err /= 0 ) then
          call handleError("do_preforce()", "Error scattering rho_col into rho")          
       endif
       rho(1:nlocal) = rho(1:nlocal) + rho_tmp(1:nlocal)
#endif
       

       do atom = 1, nlocal
          c_ident1 = c_idents_local(atom)
          
          
          call doPreforceInteraction(c_ident1, rho(atom), frho(atom), dfrhodrho(atom)) 
          pot(METALLIC_POT) = pot(METALLIC_POT) + frho(atom)
          particle_pot(atom) = particle_pot(atom) + frho(atom)
       end do

#ifdef IS_MPI
    !! communicate f(rho) and derivatives back into row and column arrays
       call gather(frho,frho_row,plan_atom_row, sc_err)
       if (sc_err /=  0) then
          call handleError("do_preforce()","MPI gather frho_row failure")
       endif
       call gather(dfrhodrho,dfrhodrho_row,plan_atom_row, sc_err)
       if (sc_err /=  0) then
          call handleError("do_preforce()","MPI gather dfrhodrho_row failure")
       endif
       call gather(frho,frho_col,plan_atom_col, sc_err)
       if (sc_err /=  0) then
          call handleError("do_preforce()","MPI gather frho_col failure")
       endif
       call gather(dfrhodrho,dfrhodrho_col,plan_atom_col, sc_err)
       if (sc_err /=  0) then
          call handleError("do_preforce()","MPI gather dfrhodrho_col failure")
       endif
#endif
       
    end if
  end subroutine do_preforce


  subroutine get_interatomic_vector(q_i, q_j, d, r_sq)

    real (kind = dp), dimension(3) :: q_i
    real (kind = dp), dimension(3) :: q_j
    real ( kind = dp ), intent(out) :: r_sq
    real( kind = dp ) :: d(3), scaled(3)
    integer i

    d(1) = q_j(1) - q_i(1)
    d(2) = q_j(2) - q_i(2)
    d(3) = q_j(3) - q_i(3)

    ! Wrap back into periodic box if necessary
    if ( SIM_uses_PBC ) then

       if( .not.boxIsOrthorhombic ) then
          ! calc the scaled coordinates.
          ! scaled = matmul(HmatInv, d)

          scaled(1) = HmatInv(1,1)*d(1) + HmatInv(1,2)*d(2) + HmatInv(1,3)*d(3)
          scaled(2) = HmatInv(2,1)*d(1) + HmatInv(2,2)*d(2) + HmatInv(2,3)*d(3)
          scaled(3) = HmatInv(3,1)*d(1) + HmatInv(3,2)*d(2) + HmatInv(3,3)*d(3)
          
          ! wrap the scaled coordinates

          scaled(1) = scaled(1) - anint(scaled(1), kind=dp)
          scaled(2) = scaled(2) - anint(scaled(2), kind=dp)
          scaled(3) = scaled(3) - anint(scaled(3), kind=dp)

          ! calc the wrapped real coordinates from the wrapped scaled 
          ! coordinates
          ! d = matmul(Hmat,scaled)
          d(1)= Hmat(1,1)*scaled(1) + Hmat(1,2)*scaled(2) + Hmat(1,3)*scaled(3)
          d(2)= Hmat(2,1)*scaled(1) + Hmat(2,2)*scaled(2) + Hmat(2,3)*scaled(3)
          d(3)= Hmat(3,1)*scaled(1) + Hmat(3,2)*scaled(2) + Hmat(3,3)*scaled(3)

       else
          ! calc the scaled coordinates.

          scaled(1) = d(1) * HmatInv(1,1)
          scaled(2) = d(2) * HmatInv(2,2)
          scaled(3) = d(3) * HmatInv(3,3)
          
          ! wrap the scaled coordinates
          
          scaled(1) = scaled(1) - anint(scaled(1), kind=dp)
          scaled(2) = scaled(2) - anint(scaled(2), kind=dp)
          scaled(3) = scaled(3) - anint(scaled(3), kind=dp)

          ! calc the wrapped real coordinates from the wrapped scaled 
          ! coordinates

          d(1) = scaled(1)*Hmat(1,1)
          d(2) = scaled(2)*Hmat(2,2)
          d(3) = scaled(3)*Hmat(3,3)

       endif

    endif

    r_sq = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)

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
    ppot_Temp = 0.0_dp

    frho_row = 0.0_dp
    frho_col = 0.0_dp
    rho_row  = 0.0_dp
    rho_col  = 0.0_dp
    rho_tmp  = 0.0_dp
    dfrhodrho_row = 0.0_dp
    dfrhodrho_col = 0.0_dp
    
#endif
    rho = 0.0_dp
    frho = 0.0_dp
    dfrhodrho = 0.0_dp

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
    unique_id_2 = AtomColToGlobal(atom2)
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
#else
    !! in the normal loop, the atom numbers are unique
    unique_id_1 = atom1
    unique_id_2 = atom2
#endif

#ifdef IS_MPI    
    do i = 1, nSkipsForRowAtom(atom1)
       if (skipsForRowAtom(atom1, i) .eq. unique_id_2) then
          skip_it = .true.
          return
       endif
    end do
#else
    do i = 1, nSkipsForLocalAtom(atom1)
       if (skipsForLocalAtom(atom1, i) .eq. unique_id_2) then
          skip_it = .true.
          return
       endif
    end do
#endif

    return
  end function skipThisPair

  function getTopoDistance(atom1, atom2) result(topoDist)
    integer, intent(in) :: atom1
    integer, intent(in) :: atom2
    integer :: topoDist
    integer :: unique_id_2
    integer :: i

#ifdef IS_MPI
    unique_id_2 = AtomColToGlobal(atom2)
#else
    unique_id_2 = atom2
#endif

    ! zero is default for unconnected (i.e. normal) pair interactions

    topoDist = 0

    do i = 1, nTopoPairsForAtom(atom1)
       if (toposForAtom(atom1, i) .eq. unique_id_2) then
          topoDist = topoDistance(atom1, i)
          return
       endif
    end do

    return
  end function getTopoDistance

  function FF_UsesDirectionalAtoms() result(doesit)
    logical :: doesit
    doesit = FF_uses_DirectionalAtoms
  end function FF_UsesDirectionalAtoms

  function FF_RequiresPrepairCalc() result(doesit)
    logical :: doesit
    doesit = FF_uses_EAM .or. FF_uses_SC
  end function FF_RequiresPrepairCalc

#ifdef PROFILE
  function getforcetime() result(totalforcetime)
    real(kind=dp) :: totalforcetime
    totalforcetime = forcetime
  end function getforcetime
#endif

  !! This cleans componets of force arrays belonging only to fortran

  subroutine add_stress_tensor(dpair, fpair, tau)

    real( kind = dp ), dimension(3), intent(in) :: dpair, fpair
    real( kind = dp ), dimension(9), intent(inout) :: tau

    ! because the d vector is the rj - ri vector, and
    ! because fx, fy, fz are the force on atom i, we need a
    ! negative sign here:  

    tau(1) = tau(1) - dpair(1) * fpair(1)
    tau(2) = tau(2) - dpair(1) * fpair(2)
    tau(3) = tau(3) - dpair(1) * fpair(3)
    tau(4) = tau(4) - dpair(2) * fpair(1)
    tau(5) = tau(5) - dpair(2) * fpair(2)
    tau(6) = tau(6) - dpair(2) * fpair(3)
    tau(7) = tau(7) - dpair(3) * fpair(1)
    tau(8) = tau(8) - dpair(3) * fpair(2)
    tau(9) = tau(9) - dpair(3) * fpair(3)

  end subroutine add_stress_tensor

end module doForces
