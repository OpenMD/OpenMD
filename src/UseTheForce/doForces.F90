!! doForces.F90
!! module doForces
!! Calculates Long Range forces.

!! @author Charles F. Vardeman II
!! @author Matthew Meineke
!! @version $Id: doForces.F90,v 1.3 2004-10-22 21:20:53 gezelter Exp $, $Date: 2004-10-22 21:20:53 $, $Name: not supported by cvs2svn $, $Revision: 1.3 $

module doForces
  use force_globals
  use simulation
  use definitions
  use atype_module
  use switcheroo
  use neighborLists  
  use lj
  use sticky_pair
  use dipole_dipole
  use charge_charge
  use reaction_field
  use gb_pair
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

  INTEGER, PARAMETER:: PREPAIR_LOOP = 1
  INTEGER, PARAMETER:: PAIR_LOOP    = 2

  logical, save :: haveRlist = .false.
  logical, save :: haveNeighborList = .false.
  logical, save :: haveSIMvariables = .false.
  logical, save :: havePropertyMap = .false.
  logical, save :: haveSaneForceField = .false.
  
  logical, save :: FF_uses_DirectionalAtoms
  logical, save :: FF_uses_LennardJones
  logical, save :: FF_uses_Electrostatic
  logical, save :: FF_uses_charges
  logical, save :: FF_uses_dipoles
  logical, save :: FF_uses_sticky
  logical, save :: FF_uses_GayBerne
  logical, save :: FF_uses_EAM
  logical, save :: FF_uses_Shapes
  logical, save :: FF_uses_FLARB
  logical, save :: FF_uses_RF

  logical, save :: SIM_uses_DirectionalAtoms
  logical, save :: SIM_uses_LennardJones
  logical, save :: SIM_uses_Electrostatics
  logical, save :: SIM_uses_Charges
  logical, save :: SIM_uses_Dipoles
  logical, save :: SIM_uses_Sticky
  logical, save :: SIM_uses_GayBerne
  logical, save :: SIM_uses_EAM
  logical, save :: SIM_uses_Shapes
  logical, save :: SIM_uses_FLARB
  logical, save :: SIM_uses_RF
  logical, save :: SIM_requires_postpair_calc
  logical, save :: SIM_requires_prepair_calc
  logical, save :: SIM_uses_PBC
  logical, save :: SIM_uses_molecular_cutoffs

  real(kind=dp), save :: rlist, rlistsq

  public :: init_FF
  public :: do_force_loop
  public :: setRlistDF

#ifdef PROFILE
  public :: getforcetime
  real, save :: forceTime = 0
  real :: forceTimeInitial, forceTimeFinal
  integer :: nLoops
#endif

  type :: Properties
     logical :: is_Directional   = .false.
     logical :: is_LennardJones  = .false.
     logical :: is_Electrostatic = .false.
     logical :: is_Charge        = .false.
     logical :: is_Dipole        = .false.
     logical :: is_Sticky        = .false.
     logical :: is_GayBerne      = .false.
     logical :: is_EAM           = .false.
     logical :: is_Shape         = .false.
     logical :: is_FLARB         = .false.
  end type Properties

  type(Properties), dimension(:),allocatable :: PropertyMap

contains

  subroutine setRlistDF( this_rlist )
    
    real(kind=dp) :: this_rlist

    rlist = this_rlist
    rlistsq = rlist * rlist
    
    haveRlist = .true.

  end subroutine setRlistDF    

  subroutine createPropertyMap(status)
    integer :: nAtypes
    integer :: status
    integer :: i
    logical :: thisProperty
    real (kind=DP) :: thisDPproperty

    status = 0

    nAtypes = getSize(atypes)

    if (nAtypes == 0) then
       status = -1
       return
    end if
        
    if (.not. allocated(PropertyMap)) then
       allocate(PropertyMap(nAtypes))
    endif

    do i = 1, nAtypes
       call getElementProperty(atypes, i, "is_Directional", thisProperty)
       PropertyMap(i)%is_Directional = thisProperty

       call getElementProperty(atypes, i, "is_LennardJones", thisProperty)
       PropertyMap(i)%is_LennardJones = thisProperty
       
       call getElementProperty(atypes, i, "is_Electrostatic", thisProperty)
       PropertyMap(i)%is_Electrostatic = thisProperty

       call getElementProperty(atypes, i, "is_Charge", thisProperty)
       PropertyMap(i)%is_Charge = thisProperty
       
       call getElementProperty(atypes, i, "is_Dipole", thisProperty)
       PropertyMap(i)%is_Dipole = thisProperty

       call getElementProperty(atypes, i, "is_Sticky", thisProperty)
       PropertyMap(i)%is_Sticky = thisProperty

       call getElementProperty(atypes, i, "is_GayBerne", thisProperty)
       PropertyMap(i)%is_GayBerne = thisProperty

       call getElementProperty(atypes, i, "is_EAM", thisProperty)
       PropertyMap(i)%is_EAM = thisProperty

       call getElementProperty(atypes, i, "is_Shape", thisProperty)
       PropertyMap(i)%is_Shape = thisProperty

       call getElementProperty(atypes, i, "is_FLARB", thisProperty)
       PropertyMap(i)%is_FLARB = thisProperty
    end do

    havePropertyMap = .true.

  end subroutine createPropertyMap

  subroutine setSimVariables()
    SIM_uses_DirectionalAtoms = SimUsesDirectionalAtoms()
    SIM_uses_LennardJones = SimUsesLennardJones()
    SIM_uses_Electrostatics = SimUsesElectrostatics()
    SIM_uses_Charges = SimUsesCharges()
    SIM_uses_Dipoles = SimUsesDipoles()
    SIM_uses_Sticky = SimUsesSticky()
    SIM_uses_GayBerne = SimUsesGayBerne()
    SIM_uses_EAM = SimUsesEAM()
    SIM_uses_Shapes = SimUsesShapes()
    SIM_uses_FLARB = SimUsesFLARB()
    SIM_uses_RF = SimUsesRF()
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
    
    if (.not. havePropertyMap) then 

       myStatus = 0

       call createPropertyMap(myStatus)

       if (myStatus .ne. 0) then
          write(default_error, *) 'createPropertyMap failed in doForces!'
          error = -1
          return
       endif
    endif

    if (.not. haveSIMvariables) then
       call setSimVariables()
    endif

    if (.not. haveRlist) then
       write(default_error, *) 'rList has not been set in doForces!'
       error = -1
       return
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
    

  subroutine init_FF(use_RF_c, thisStat)

    logical, intent(in) :: use_RF_c

    integer, intent(out) :: thisStat   
    integer :: my_status, nMatches
    integer, pointer :: MatchList(:) => null()
    real(kind=dp) :: rcut, rrf, rt, dielect

    !! assume things are copacetic, unless they aren't
    thisStat = 0

    !! Fortran's version of a cast:
    FF_uses_RF = use_RF_c
    
    !! init_FF is called *after* all of the atom types have been 
    !! defined in atype_module using the new_atype subroutine.
    !!
    !! this will scan through the known atypes and figure out what
    !! interactions are used by the force field.     
  
    FF_uses_DirectionalAtoms = .false.
    FF_uses_LennardJones = .false.
    FF_uses_Electrostatic = .false.
    FF_uses_Charges = .false.    
    FF_uses_Dipoles = .false.
    FF_uses_Sticky = .false.
    FF_uses_GayBerne = .false.
    FF_uses_EAM = .false.
    FF_uses_Shapes = .false.
    FF_uses_FLARB = .false.
    
    call getMatchingElementList(atypes, "is_Directional", .true., &
         nMatches, MatchList)
    if (nMatches .gt. 0) FF_uses_DirectionalAtoms = .true.

    call getMatchingElementList(atypes, "is_LennardJones", .true., &
         nMatches, MatchList)
    if (nMatches .gt. 0) FF_uses_LennardJones = .true.
    
    call getMatchingElementList(atypes, "is_Electrostatic", .true., &
         nMatches, MatchList)
    if (nMatches .gt. 0) then
       FF_uses_Electrostatic = .true.
    endif

    call getMatchingElementList(atypes, "is_Charge", .true., &
         nMatches, MatchList)
    if (nMatches .gt. 0) then
       FF_uses_charges = .true.   
       FF_uses_electrostatic = .true.
    endif
    
    call getMatchingElementList(atypes, "is_Dipole", .true., &
         nMatches, MatchList)
    if (nMatches .gt. 0) then 
       FF_uses_dipoles = .true.
       FF_uses_electrostatic = .true.
       FF_uses_DirectionalAtoms = .true.
    endif
    
    call getMatchingElementList(atypes, "is_Sticky", .true., nMatches, &
         MatchList) 
    if (nMatches .gt. 0) then
       FF_uses_Sticky = .true.
       FF_uses_DirectionalAtoms = .true.
    endif
    
    call getMatchingElementList(atypes, "is_GayBerne", .true., &
         nMatches, MatchList)
    if (nMatches .gt. 0) then
       FF_uses_GayBerne = .true.
       FF_uses_DirectionalAtoms = .true.
    endif
    
    call getMatchingElementList(atypes, "is_EAM", .true., nMatches, MatchList)
    if (nMatches .gt. 0) FF_uses_EAM = .true.
    
    call getMatchingElementList(atypes, "is_Shape", .true., &
         nMatches, MatchList)
    if (nMatches .gt. 0) then
       FF_uses_Shapes = .true.
       FF_uses_DirectionalAtoms = .true.
    endif

    call getMatchingElementList(atypes, "is_FLARB", .true., &
         nMatches, MatchList)
    if (nMatches .gt. 0) FF_uses_FLARB = .true.

    !! Assume sanity (for the sake of argument)
    haveSaneForceField = .true.
    
    !! check to make sure the FF_uses_RF setting makes sense
    
    if (FF_uses_dipoles) then
       if (FF_uses_RF) then
          dielect = getDielect()
          call initialize_rf(dielect)
       endif
    else
       if (FF_uses_RF) then          
          write(default_error,*) 'Using Reaction Field with no dipoles?  Huh?'
          thisStat = -1
          haveSaneForceField = .false.
          return
       endif
    endif 

    if (FF_uses_sticky) then
       call check_sticky_FF(my_status)
       if (my_status /= 0) then
          thisStat = -1
          haveSaneForceField = .false.
          return
       end if
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

    if (FF_uses_GayBerne) then
       call check_gb_pair_FF(my_status)
       if (my_status .ne. 0) then
          thisStat = -1
          haveSaneForceField = .false.
          return
       endif
    endif

    if (FF_uses_GayBerne .and. FF_uses_LennardJones) then
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
  subroutine do_force_loop(q, q_group, A, u_l, f, t, tau, pot, &
       do_pot_c, do_stress_c, error)
    !! Position array provided by C, dimensioned by getNlocal
    real ( kind = dp ), dimension(3, nLocal) :: q
    !! molecular center-of-mass position array
    real ( kind = dp ), dimension(3, nGroups) :: q_group
    !! Rotation Matrix for each long range particle in simulation.
    real( kind = dp), dimension(9, nLocal) :: A    
    !! Unit vectors for dipoles (lab frame)
    real( kind = dp ), dimension(3,nLocal) :: u_l
    !! Force array provided by C, dimensioned by getNlocal
    real ( kind = dp ), dimension(3,nLocal) :: f
    !! Torsion array provided by C, dimensioned by getNlocal
    real( kind = dp ), dimension(3,nLocal) :: t    

    !! Stress Tensor
    real( kind = dp), dimension(9) :: tau   
    real ( kind = dp ) :: pot
    logical ( kind = 2) :: do_pot_c, do_stress_c
    logical :: do_pot
    logical :: do_stress
    logical :: in_switching_region
#ifdef IS_MPI 
    real( kind = DP ) :: pot_local
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

    real(kind=dp) :: listSkin = 1.0   
    
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
       call gather(u_l, u_l_Row, plan_atom_row_3d)
       call gather(u_l, u_l_Col, plan_atom_col_3d)
       
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
             call get_interatomic_vector(q_group_Row(:,i), &
                  q_group_Col(:,j), d_grp, rgrpsq)
#else
             call get_interatomic_vector(q_group(:,i), &
                  q_group(:,j), d_grp, rgrpsq)
#endif

             if (rgrpsq < rlistsq) then
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
                              u_l, A, f, t, pot_local)
#else
                         call do_prepair(atom1, atom2, ratmsq, d_atm, sw, &
                              rgrpsq, d_grp, do_pot, do_stress, &
                              u_l, A, f, t, pot)
#endif                                               
                      else
#ifdef IS_MPI                      
                         call do_pair(atom1, atom2, ratmsq, d_atm, sw, &
                              do_pot, &
                              u_l, A, f, t, pot_local, vpair, fpair)
#else
                         call do_pair(atom1, atom2, ratmsq, d_atm, sw, &
                              do_pot,  &
                              u_l, A, f, t, pot, vpair, fpair)
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
       call scatter(pot_Row, pot_Temp, plan_atom_row)
       
       ! scatter/gather pot_local into all other procs
       ! add resultant to get total pot
       do i = 1, nlocal
          pot_local = pot_local + pot_Temp(i)
       enddo
       
       pot_Temp = 0.0_DP 
       
       call scatter(pot_Col, pot_Temp, plan_atom_col)
       do i = 1, nlocal
          pot_local = pot_local + pot_Temp(i)
       enddo
       
    endif
#endif
    
    if (FF_RequiresPostpairCalc() .and. SIM_requires_postpair_calc) then
       
       if (FF_uses_RF .and. SIM_uses_RF) then
          
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
             
             if (PropertyMap(me_i)%is_Dipole) then
                
                mu_i = getDipoleMoment(me_i)
                
                !! The reaction field needs to include a self contribution 
                !! to the field:
                call accumulate_self_rf(i, mu_i, u_l)
                !! Get the reaction field contribution to the 
                !! potential and torques:
                call reaction_field_final(i, mu_i, u_l, rfpot, t, do_pot)
#ifdef IS_MPI
                pot_local = pot_local + rfpot
#else
                pot = pot + rfpot
      
#endif
             endif             
          enddo
       endif
    endif
    
    
#ifdef IS_MPI
    
    if (do_pot) then
       pot = pot + pot_local
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
       u_l, A, f, t, pot, vpair, fpair)

    real( kind = dp ) :: pot, vpair, sw
    real( kind = dp ), dimension(3) :: fpair
    real( kind = dp ), dimension(nLocal)   :: mfact
    real( kind = dp ), dimension(3,nLocal) :: u_l
    real( kind = dp ), dimension(9,nLocal) :: A
    real( kind = dp ), dimension(3,nLocal) :: f
    real( kind = dp ), dimension(3,nLocal) :: t

    logical, intent(inout) :: do_pot
    integer, intent(in) :: i, j
    real ( kind = dp ), intent(inout) :: rijsq
    real ( kind = dp )                :: r
    real ( kind = dp ), intent(inout) :: d(3)
    integer :: me_i, me_j

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
    
    if (FF_uses_LennardJones .and. SIM_uses_LennardJones) then
       
       if ( PropertyMap(me_i)%is_LennardJones .and. &
            PropertyMap(me_j)%is_LennardJones ) then
          call do_lj_pair(i, j, d, r, rijsq, sw, vpair, fpair, pot, f, do_pot)
       endif
       
    endif
    
    if (FF_uses_charges .and. SIM_uses_charges) then
       
       if (PropertyMap(me_i)%is_Charge .and. PropertyMap(me_j)%is_Charge) then
          call do_charge_pair(i, j, d, r, rijsq, sw, vpair, fpair, &
               pot, f, do_pot)
       endif
       
    endif
    
    if (FF_uses_dipoles .and. SIM_uses_dipoles) then
       
       if ( PropertyMap(me_i)%is_Dipole .and. PropertyMap(me_j)%is_Dipole) then
          call do_dipole_pair(i, j, d, r, rijsq, sw, vpair, fpair, &
               pot, u_l, f, t, do_pot)
          if (FF_uses_RF .and. SIM_uses_RF) then
             call accumulate_rf(i, j, r, u_l, sw)
             call rf_correct_forces(i, j, d, r, u_l, sw, f, fpair)
          endif
       endif

    endif

    if (FF_uses_Sticky .and. SIM_uses_sticky) then

       if ( PropertyMap(me_i)%is_Sticky .and. PropertyMap(me_j)%is_Sticky) then
          call do_sticky_pair(i, j, d, r, rijsq, sw, vpair, fpair, &
               pot, A, f, t, do_pot)
       endif
       
    endif


    if (FF_uses_GayBerne .and. SIM_uses_GayBerne) then
       
       if ( PropertyMap(me_i)%is_GayBerne .and. &
            PropertyMap(me_j)%is_GayBerne) then
          call do_gb_pair(i, j, d, r, rijsq, sw, vpair, fpair, &
               pot, u_l, f, t, do_pot)
       endif
       
    endif
    
    if (FF_uses_EAM .and. SIM_uses_EAM) then
       
       if ( PropertyMap(me_i)%is_EAM .and. PropertyMap(me_j)%is_EAM) then
          call do_eam_pair(i, j, d, r, rijsq, sw, vpair, fpair, pot, f, &
               do_pot)
       endif
       
    endif

    if (FF_uses_Shapes .and. SIM_uses_Shapes) then
       
       if ( PropertyMap(me_i)%is_Shape .and. &
            PropertyMap(me_j)%is_Shape ) then
          call do_shape_pair(i, j, d, r, rijsq, sw, vpair, fpair, &
               pot, u_l, f, t, do_pot)
       endif
       
    endif
    
  end subroutine do_pair

  subroutine do_prepair(i, j, rijsq, d, sw, rcijsq, dc, &
       do_pot, do_stress, u_l, A, f, t, pot)

   real( kind = dp ) :: pot, sw
   real( kind = dp ), dimension(3,nLocal) :: u_l
   real (kind=dp), dimension(9,nLocal) :: A
   real (kind=dp), dimension(3,nLocal) :: f
   real (kind=dp), dimension(3,nLocal) :: t
   
   logical, intent(inout) :: do_pot, do_stress
   integer, intent(in) :: i, j
   real ( kind = dp ), intent(inout)    :: rijsq, rcijsq
   real ( kind = dp )                :: r, rc
   real ( kind = dp ), intent(inout) :: d(3), dc(3)
   
   logical :: is_EAM_i, is_EAM_j
   
   integer :: me_i, me_j
   

    r = sqrt(rijsq)
    if (SIM_uses_molecular_cutoffs) then
       rc = sqrt(rcijsq)
    else
       rc = r
    endif
   

#ifdef IS_MPI   
   me_i = atid_row(i)
   me_j = atid_col(j)   
#else   
   me_i = atid(i)
   me_j = atid(j)   
#endif
   
   if (FF_uses_EAM .and. SIM_uses_EAM) then
      
      if (PropertyMap(me_i)%is_EAM .and. PropertyMap(me_j)%is_EAM) &
           call calc_EAM_prepair_rho(i, j, d, r, rijsq )
      
   endif
   
 end subroutine do_prepair
 
 
 subroutine do_preforce(nlocal,pot)
   integer :: nlocal
   real( kind = dp ) :: pot
   
   if (FF_uses_EAM .and. SIM_uses_EAM) then
      call calc_EAM_preforce_Frho(nlocal,pot)
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
   
   u_l_Row = 0.0_dp
   u_l_Col = 0.0_dp
   
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
   doesit = FF_uses_DirectionalAtoms .or. FF_uses_Dipoles .or. &
        FF_uses_Sticky .or. FF_uses_GayBerne .or. FF_uses_Shapes
 end function FF_UsesDirectionalAtoms
 
 function FF_RequiresPrepairCalc() result(doesit)
   logical :: doesit
   doesit = FF_uses_EAM
 end function FF_RequiresPrepairCalc
 
 function FF_RequiresPostpairCalc() result(doesit)
   logical :: doesit
   doesit = FF_uses_RF
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

!! Interfaces for C programs to module....

 subroutine initFortranFF(use_RF_c, thisStat)
    use doForces, ONLY: init_FF
    logical, intent(in) :: use_RF_c

    integer, intent(out) :: thisStat   
    call init_FF(use_RF_c, thisStat)

 end subroutine initFortranFF

  subroutine doForceloop(q, q_group, A, u_l, f, t, tau, pot, &
       do_pot_c, do_stress_c, error)
       
       use definitions, ONLY: dp
       use simulation
       use doForces, ONLY: do_force_loop
    !! Position array provided by C, dimensioned by getNlocal
    real ( kind = dp ), dimension(3, nLocal) :: q
    !! molecular center-of-mass position array
    real ( kind = dp ), dimension(3, nGroups) :: q_group
    !! Rotation Matrix for each long range particle in simulation.
    real( kind = dp), dimension(9, nLocal) :: A    
    !! Unit vectors for dipoles (lab frame)
    real( kind = dp ), dimension(3,nLocal) :: u_l
    !! Force array provided by C, dimensioned by getNlocal
    real ( kind = dp ), dimension(3,nLocal) :: f
    !! Torsion array provided by C, dimensioned by getNlocal
    real( kind = dp ), dimension(3,nLocal) :: t    

    !! Stress Tensor
    real( kind = dp), dimension(9) :: tau   
    real ( kind = dp ) :: pot
    logical ( kind = 2) :: do_pot_c, do_stress_c
    integer :: error
    
    call do_force_loop(q, q_group, A, u_l, f, t, tau, pot, &
       do_pot_c, do_stress_c, error)
       
 end subroutine doForceloop
