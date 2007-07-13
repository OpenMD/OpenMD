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
!! @version $Id: MetalNonMetal.F90,v 1.1 2007-07-13 14:21:07 chuckv Exp $, $Date: 2007-07-13 14:21:07 $, $Name: not supported by cvs2svn $, $Revision: 1.1 $


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


contains



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
       write(*,*) 'unknown MnM interaction type!'
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

end module MetalNonMetal
