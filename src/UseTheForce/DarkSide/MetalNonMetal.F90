!!
!! Copyright (c) 2007 The University of Notre Dame. All Rights Reserved.
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


!! Calculates Metal-Non Metal interactions.
!! @author Charles F. Vardeman II 
!! @version $Id: MetalNonMetal.F90,v 1.2 2007-07-17 18:54:47 chuckv Exp $, $Date: 2007-07-17 18:54:47 $, $Name: not supported by cvs2svn $, $Revision: 1.2 $


module MetalNonMetal
  use definitions
  use atype_module
  use vector_class
  use simulation
  use status
  use fForceOptions
#ifdef IS_MPI
  use mpiSimulation
#endif
  use force_globals

  implicit none
  PRIVATE
#define __FORTRAN90
#include "UseTheForce/DarkSide/fInteractionMap.h"
#include "UseTheForce/DarkSide/fMnMInteractions.h"

  logical, save :: useGeometricDistanceMixing = .false.
  logical, save :: haveInteractionLookup = .false.

  real(kind=DP), save :: defaultCutoff = 0.0_DP
  logical, save :: defaultShiftPot = .false.
  logical, save :: defaultShiftFrc = .false.
  logical, save :: haveDefaultCutoff = .false.

  type :: MnMinteraction
     integer :: metal_atid
     integer :: nonmetal_atid
     integer :: interaction_type
     real(kind=dp) :: R0
     real(kind=dp) :: D0
     real(kind=dp) :: beta0
     real(kind=dp) :: betaH
     real(kind=dp) :: alpha
     real(kind=dp) :: gamma     
     real(kind=dp) :: sigma
     real(kind=dp) :: epsilon
     real(kind=dp) :: rCut = 0.0_dp
     logical       :: rCutWasSet = .false.
     logical       :: shiftedPot = .false.
     logical       :: shiftedFrc = .false.
  end type MnMinteraction

  type :: MnMinteractionMap
     PRIVATE
     integer :: initialCapacity = 10
     integer :: capacityIncrement = 0
     integer :: interactionCount = 0
     type(MnMinteraction), pointer :: interactions(:) => null()
  end type MnMinteractionMap

  type (MnMInteractionMap), pointer :: MnM_Map

  integer,  allocatable, dimension(:,:) :: MnMinteractionLookup

  public :: setMnMDefaultCutoff
  public :: addInteraction
  public :: deleteInteractions
  public :: MNMtype
  public :: do_mnm_pair

contains


  subroutine do_mnm_pair(Atom1, Atom2, D, Rij, R2, Rcut, Sw, Vpair, Fpair, &
       Pot, A, F,t, Do_pot)
    integer, intent(inout) ::  atom1, atom2
    integer :: atid1, atid2, ljt1, ljt2
    real( kind = dp ), intent(inout) :: rij, r2, rcut
    real( kind = dp ), intent(inout) :: pot, sw, vpair
    real( kind = dp ), intent(inout), dimension(3,nLocal) :: f 
    real (kind=dp),intent(inout), dimension(9,nLocal) :: A
    real (kind=dp),intent(inout), dimension(3,nLocal) :: t   
    real( kind = dp ), intent(inout), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair
    logical, intent(inout) :: do_pot

    integer :: interaction_id
    integer :: interaction_type

#ifdef IS_MPI
    atid1 = atid_Row(atom1)
    atid2 = atid_Col(atom2)
#else
    atid1 = atid(atom1)
    atid2 = atid(atom2)
#endif

    interaction_id = MnMinteractionLookup(atid1, atid2)
    interaction_type = MnM_Map%interactions(interaction_id)%interaction_type

    select case (interaction_type)
    case (MNM_LENNARDJONES)
       call calc_mnm_lennardjones(Atom1, Atom2, D, Rij, R2, Rcut, Sw, Vpair, Fpair, &
       Pot, F, Do_pot, interaction_id)
    case(MNM_REPULSIVEMORSE, MNM_SHIFTEDMORSE)
       call calc_mnm_morse(Atom1, Atom2, D, Rij, R2, Rcut, Sw, Vpair, Fpair, &
       Pot, F, Do_pot, interaction_id,interaction_type)
    case(MNM_MAW)
       call calc_mnm_maw(Atom1, Atom2, D, Rij, R2, Rcut, Sw, Vpair, Fpair, &
       Pot,A, F,t, Do_pot, interaction_id)
    end select

  end subroutine do_mnm_pair

  subroutine calc_mnm_lennardjones(Atom1, Atom2, D, Rij, R2, Rcut, Sw, Vpair, Fpair, &
       Pot, F, Do_pot, interaction_id)
    
    integer, intent(inout) ::  atom1, atom2
    real( kind = dp ), intent(inout) :: rij, r2, rcut
    real( kind = dp ), intent(inout) :: pot, sw, vpair
    real( kind = dp ), intent(inout), dimension(3,nLocal) :: f    
    real( kind = dp ), intent(inout), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair
    logical, intent(inout) :: do_pot
    integer, intent(in) :: interaction_id

    ! local Variables
    real( kind = dp ) :: drdx, drdy, drdz
    real( kind = dp ) :: fx, fy, fz
    real( kind = dp ) :: myPot, myPotC, myDeriv, myDerivC, ros, rcos
    real( kind = dp ) :: pot_temp, dudr
    real( kind = dp ) :: sigmai
    real( kind = dp ) :: epsilon
    logical :: isSoftCore, shiftedPot, shiftedFrc
    integer :: id1, id2, localError


    sigmai     = MnM_Map%interactions(interaction_id)%sigma
    epsilon    = MnM_Map%interactions(interaction_id)%epsilon
    shiftedPot = MnM_Map%interactions(interaction_id)%shiftedPot
    shiftedFrc = MnM_Map%interactions(interaction_id)%shiftedFrc

    ros = rij * sigmai


    call getLJfunc(ros, myPot, myDeriv)
    
    if (shiftedPot) then
       rcos = rcut * sigmai
       call getLJfunc(rcos, myPotC, myDerivC) 
       myDerivC = 0.0_dp
    elseif (shiftedFrc) then
       rcos = rcut * sigmai
       call getLJfunc(rcos, myPotC, myDerivC)
       myPotC = myPotC + myDerivC * (rij - rcut) * sigmai
    else
       myPotC = 0.0_dp
       myDerivC = 0.0_dp
    endif
    
 

    pot_temp = epsilon * (myPot - myPotC)
    vpair = vpair + pot_temp
    dudr = sw * epsilon * (myDeriv - myDerivC) * sigmai

    drdx = d(1) / rij
    drdy = d(2) / rij
    drdz = d(3) / rij

    fx = dudr * drdx
    fy = dudr * drdy
    fz = dudr * drdz

#ifdef IS_MPI
    if (do_pot) then
       pot_Row(VDW_POT,atom1) = pot_Row(VDW_POT,atom1) + sw*pot_temp*0.5
       pot_Col(VDW_POT,atom2) = pot_Col(VDW_POT,atom2) + sw*pot_temp*0.5
    endif

    f_Row(1,atom1) = f_Row(1,atom1) + fx 
    f_Row(2,atom1) = f_Row(2,atom1) + fy
    f_Row(3,atom1) = f_Row(3,atom1) + fz

    f_Col(1,atom2) = f_Col(1,atom2) - fx 
    f_Col(2,atom2) = f_Col(2,atom2) - fy
    f_Col(3,atom2) = f_Col(3,atom2) - fz       

#else
    if (do_pot) pot = pot + sw*pot_temp

    f(1,atom1) = f(1,atom1) + fx
    f(2,atom1) = f(2,atom1) + fy
    f(3,atom1) = f(3,atom1) + fz

    f(1,atom2) = f(1,atom2) - fx
    f(2,atom2) = f(2,atom2) - fy
    f(3,atom2) = f(3,atom2) - fz
#endif

#ifdef IS_MPI
    id1 = AtomRowToGlobal(atom1)
    id2 = AtomColToGlobal(atom2)
#else
    id1 = atom1
    id2 = atom2
#endif

    if (molMembershipList(id1) .ne. molMembershipList(id2)) then

       fpair(1) = fpair(1) + fx
       fpair(2) = fpair(2) + fy
       fpair(3) = fpair(3) + fz

    endif

    return    


  end subroutine calc_mnm_lennardjones






  subroutine calc_mnm_morse(Atom1, Atom2, D, Rij, R2, Rcut, Sw, Vpair, Fpair, &
       Pot, f, Do_pot, interaction_id, interaction_type)
    integer, intent(inout) ::  atom1, atom2
    real( kind = dp ), intent(inout) :: rij, r2, rcut
    real( kind = dp ), intent(inout) :: pot, sw, vpair
    real( kind = dp ), intent(inout), dimension(3,nLocal) :: f    
    real( kind = dp ), intent(inout), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair
    logical, intent(inout) :: do_pot
    integer, intent(in) :: interaction_id, interaction_type

    ! Local Variables
    real(kind=dp) :: Beta0
    real(kind=dp) :: R0
    real(kind=dp) :: D0
    real(kind=dp) :: expt
    real(kind=dp) :: expt2
    real(kind=dp) :: expfnc
    real(kind=dp) :: expfnc2
    real(kind=dp) :: D_expt
    real(kind=dp) :: D_expt2
    real(kind=dp) :: rcos
    real(kind=dp) :: exptC
    real(kind=dp) :: expt2C
    real(kind=dp) :: expfncC
    real(kind=dp) :: expfnc2C
    real(kind=dp) :: D_expt2C
    real(kind=dp) :: D_exptC
    
    real(kind=dp) :: myPot
    real(kind=dp) :: myPotC
    real(kind=dp) :: myDeriv
    real(kind=dp) :: myDerivC
    real(kind=dp) :: pot_temp
    real(kind=dp) :: fx,fy,fz
    real(kind=dp) :: drdx,drdy,drdz
    real(kind=dp) :: dudr
    integer :: id1,id2


    D0     = MnM_Map%interactions(interaction_id)%D0
    R0    = MnM_Map%interactions(interaction_id)%r0
    Beta0 = MnM_Map%interactions(interaction_id)%Beta0
!    shiftedFrc = MnM_Map%interactions(interaction_id)%shiftedFrc


! V(r) = D_e exp(-a(r-re)(exp(-a(r-re))-2)

    expt     = -R0*(rij-R0) 
    expt2    = 2.0_dp*expt
    expfnc   = exp(expt)
    expfnc2  = exp(expt2)
    D_expt   = D0*expt
    D_expt2  = D0*expt2

    rcos      = rcut*Beta0
    exptC     = -R0*(rcos-R0) 
    expt2C    = 2.0_dp*expt
    expfncC   = exp(expt)
    expfnc2C  = exp(expt2)
    D_exptC   = D0*expt
    D_expt2C  = D0*expt2

    select case (interaction_type)


    case (MNM_SHIFTEDMORSE)
       
       myPot  = D_expt  * (D_expt  - 2.0_dp)
       myPotC = D_exptC * (D_exptC - 2.0_dp)

       myDeriv   = -(D_expt2  - D_expt)  * (-2.0_dp + D_expt)
       myDerivC  = -(D_expt2C - D_exptC) * (-2.0_dp + D_exptC)

    case (MNM_REPULSIVEMORSE)

       myPot  = D_expt2
       myPotC = D_expt2C

       myDeriv  = -2.0_dp * D_expt2
       myDerivC = -2.0_dp * D_expt2C

    end select
 
    myPotC = myPotC + myDerivC*(rij - rcut)*Beta0

    pot_temp = (myPot - myPotC)
    vpair = vpair + pot_temp
    dudr = sw * (myDeriv - myDerivC)

    drdx = d(1) / rij
    drdy = d(2) / rij
    drdz = d(3) / rij

    fx = dudr * drdx
    fy = dudr * drdy
    fz = dudr * drdz

#ifdef IS_MPI
    if (do_pot) then
       pot_Row(VDW_POT,atom1) = pot_Row(VDW_POT,atom1) + sw*pot_temp*0.5
       pot_Col(VDW_POT,atom2) = pot_Col(VDW_POT,atom2) + sw*pot_temp*0.5
    endif

    f_Row(1,atom1) = f_Row(1,atom1) + fx 
    f_Row(2,atom1) = f_Row(2,atom1) + fy
    f_Row(3,atom1) = f_Row(3,atom1) + fz

    f_Col(1,atom2) = f_Col(1,atom2) - fx 
    f_Col(2,atom2) = f_Col(2,atom2) - fy
    f_Col(3,atom2) = f_Col(3,atom2) - fz       

#else
    if (do_pot) pot = pot + sw*pot_temp

    f(1,atom1) = f(1,atom1) + fx
    f(2,atom1) = f(2,atom1) + fy
    f(3,atom1) = f(3,atom1) + fz

    f(1,atom2) = f(1,atom2) - fx
    f(2,atom2) = f(2,atom2) - fy
    f(3,atom2) = f(3,atom2) - fz
#endif

#ifdef IS_MPI
    id1 = AtomRowToGlobal(atom1)
    id2 = AtomColToGlobal(atom2)
#else
    id1 = atom1
    id2 = atom2
#endif

    if (molMembershipList(id1) .ne. molMembershipList(id2)) then

       fpair(1) = fpair(1) + fx
       fpair(2) = fpair(2) + fy
       fpair(3) = fpair(3) + fz

    endif

    return    

  end subroutine calc_mnm_morse


  

  subroutine calc_mnm_maw(Atom1, Atom2, D, Rij, R2, Rcut, Sw, Vpair, Fpair, &
       Pot,A, F,t, Do_pot, interaction_id)
    integer, intent(inout) ::  atom1, atom2
    real( kind = dp ), intent(inout) :: rij, r2, rcut
    real( kind = dp ), intent(inout) :: pot, sw, vpair
    real( kind = dp ), intent(inout), dimension(3,nLocal) :: f    
    real (kind=dp),intent(inout), dimension(9,nLocal) :: A
    real (kind=dp),intent(inout), dimension(3,nLocal) :: t   

    real( kind = dp ), intent(inout), dimension(3) :: d
    real( kind = dp ), intent(inout), dimension(3) :: fpair
    logical, intent(inout) :: do_pot

    integer, intent(in) :: interaction_id

  end subroutine calc_mnm_maw


  subroutine  setMnMDefaultCutoff(thisRcut, shiftedPot, shiftedFrc)
    real(kind=dp), intent(in) :: thisRcut
    logical, intent(in) :: shiftedPot
    logical, intent(in) :: shiftedFrc
    integer i, nInteractions
    defaultCutoff = thisRcut
    defaultShiftPot = shiftedPot
    defaultShiftFrc = shiftedFrc

    if(MnM_Map%interactionCount /= 0) then
       nInteractions = MnM_Map%interactionCount

       do i = 1, nInteractions
          MnM_Map%interactions(i)%shiftedPot = shiftedPot
          MnM_Map%interactions(i)%shiftedFrc = shiftedFrc
          MnM_Map%interactions(i)%rCut = thisRcut
          MnM_Map%interactions(i)%rCutWasSet = .true.
       enddo
    end if

  end subroutine setMnMDefaultCutoff

  subroutine copyAllData(v1, v2)
    type(MnMinteractionMap), pointer  :: v1
    type(MnMinteractionMap), pointer  :: v2
    integer :: i, j

    do i = 1, v1%interactionCount
       v2%interactions(i) = v1%interactions(i)
    enddo

    v2%interactionCount = v1%interactionCount
    return
  end subroutine copyAllData

  subroutine addInteraction(myInteraction)
    type(MNMtype) :: myInteraction
    type(MnMinteraction) :: nt
    integer :: id

    nt%interaction_type = myInteraction%MNMInteractionType
    nt%metal_atid = myInteraction%metal_atid
    nt%nonmetal_atid = myInteraction%nonmetal_atid
    
    select case (nt%interaction_type)
    case (MNM_LENNARDJONES)
       nt%sigma = myInteraction%sigma
       nt%epsilon = myInteraction%epsilon
    case(MNM_REPULSIVEMORSE, MNM_SHIFTEDMORSE)
       nt%R0 = myInteraction%R0
       nt%D0 = myInteraction%D0
       nt%beta0 = myInteraction%beta0
    case(MNM_MAW)
       nt%R0 = myInteraction%R0
       nt%D0 = myInteraction%D0
       nt%beta0 = myInteraction%beta0
       nt%betaH = myInteraction%betaH
       nt%alpha = myInteraction%alpha
       nt%gamma = myInteraction%gamma
    case default
       call handleError("MNM", "Unknown Interaction type")
    end select
    
    if (.not. associated(MnM_Map)) then
       call ensureCapacityHelper(MnM_Map, 1)
    else
       call ensureCapacityHelper(MnM_Map, MnM_Map%interactionCount + 1)
    end if
    
    MnM_Map%interactionCount = MnM_Map%interactionCount + 1
    id = MnM_Map%interactionCount
    MnM_Map%interactions(id) = nt
  end subroutine addInteraction

  subroutine ensureCapacityHelper(this, minCapacity)
    type(MnMinteractionMap), pointer :: this, that
    integer, intent(in) :: minCapacity
    integer :: oldCapacity 
    integer :: newCapacity
    logical :: resizeFlag 

    resizeFlag = .false.

    !  first time: allocate a new vector with default size

    if (.not. associated(this)) then
       this => MnMinitialize(minCapacity, 0)
    endif

    oldCapacity = size(this%interactions)

    if (minCapacity > oldCapacity) then
       if (this%capacityIncrement .gt. 0) then
          newCapacity = oldCapacity + this%capacityIncrement
       else
          newCapacity = oldCapacity * 2
       endif
       if (newCapacity .lt. minCapacity) then
          newCapacity = minCapacity
       endif
       resizeFlag = .true.
    else
       newCapacity = oldCapacity
    endif

    if (resizeFlag) then
       that => MnMinitialize(newCapacity, this%capacityIncrement)
       call copyAllData(this, that)
       this => MnMdestroy(this)
       this => that
    endif
  end subroutine ensureCapacityHelper

  function MnMinitialize(cap, capinc) result(this)
    integer, intent(in) :: cap, capinc
    integer :: error
    type(MnMinteractionMap), pointer :: this 

    nullify(this)

    if (cap < 0) then
       write(*,*) 'Bogus Capacity:', cap
       return
    endif
    allocate(this,stat=error)
    if ( error /= 0 ) then
       write(*,*) 'Could not allocate MnMinteractionMap!'
       return
    end if

    this%initialCapacity = cap
    this%capacityIncrement = capinc

    allocate(this%interactions(this%initialCapacity), stat=error)
    if(error /= 0) write(*,*) 'Could not allocate MnMinteraction!'

  end function MnMinitialize

  subroutine createInteractionLookup(this)
    type(MnMinteractionMap), pointer :: this
    integer :: biggestAtid, i, metal_atid, nonmetal_atid, error

    biggestAtid=-1
    do i = 1, this%interactionCount
       metal_atid = this%interactions(i)%metal_atid
       nonmetal_atid = this%interactions(i)%nonmetal_atid

       if (metal_atid .gt. biggestAtid) biggestAtid = metal_atid
       if (nonmetal_atid .gt. biggestAtid) biggestAtid = nonmetal_atid
    enddo

    allocate(MnMinteractionLookup(biggestAtid,biggestAtid), stat=error)
    if (error /= 0) write(*,*) 'Could not allocate MnMinteractionLookup'

    do i = 1, this%interactionCount
       metal_atid = this%interactions(i)%metal_atid
       nonmetal_atid = this%interactions(i)%nonmetal_atid
    
       MnMinteractionLookup(metal_atid, nonmetal_atid) = i
       MnMinteractionLookup(nonmetal_atid, metal_atid) = i
    enddo
  end subroutine createInteractionLookup
    

  function MnMdestroy(this) result(null_this)
    logical :: done
    type(MnMinteractionMap), pointer :: this 
    type(MnMinteractionMap), pointer :: null_this 

    if (.not. associated(this)) then
       null_this => null()
       return
    end if

    !! Walk down the list and deallocate each of the map's components
    if(associated(this%interactions)) then
       deallocate(this%interactions)
       this%interactions=>null()
    endif
    deallocate(this)
    this => null()
    null_this => null()
  end function MnMdestroy


  subroutine deleteInteractions()    
    MnM_Map => MnMdestroy(MnM_Map)
    return
  end subroutine deleteInteractions


  subroutine getLJfunc(r, myPot, myDeriv)

    real(kind=dp), intent(in) :: r
    real(kind=dp), intent(inout) :: myPot, myDeriv
    real(kind=dp) :: ri, ri2, ri6, ri7, ri12, ri13
    real(kind=dp) :: a, b, c, d, dx
    integer :: j

    ri = 1.0_DP / r
    ri2 = ri*ri
    ri6 = ri2*ri2*ri2
    ri7 = ri6*ri
    ri12 = ri6*ri6
    ri13 = ri12*ri
    
    myPot = 4.0_DP * (ri12 - ri6)
    myDeriv = 24.0_DP * (ri7 - 2.0_DP * ri13)
    
    return
  end subroutine getLJfunc

  subroutine getSoftFunc(r, myPot, myDeriv)
    
    real(kind=dp), intent(in) :: r
    real(kind=dp), intent(inout) :: myPot, myDeriv
    real(kind=dp) :: ri, ri2, ri6, ri7
    real(kind=dp) :: a, b, c, d, dx
    integer :: j
    
    ri = 1.0_DP / r    
    ri2 = ri*ri
    ri6 = ri2*ri2*ri2
    ri7 = ri6*ri
    myPot = 4.0_DP * (ri6)
    myDeriv = - 24.0_DP * ri7 
    
    return
  end subroutine getSoftFunc






end module MetalNonMetal
