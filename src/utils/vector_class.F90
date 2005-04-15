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

! vector_class.F90
!! Module Vector_class
!! Fortran 95 Vector class module. Similar to java.util vector class.
!! 
!! The Vector class implements a growable array of objects. Like an array, 
!! it contains components that can be accessed using an integer index. 
!! However, the size of a Vector can grow as needed to accommodate 
!! adding and removing items after the Vector has been created.
!! Each vector tries to optimize storage management by maintaining a 
!! capacity and a capacityIncrement. The capacity is always at least as 
!! large as the vector size; 
!! it is usually larger because as components are added to the vector, 
!! the vector's storage increases in chunks the size of capacityIncrement. 
!! An application can increase the capacity of a vector before inserting a 
!! large number of components; this reduces the amount of incremental 
!! reallocation. 
!! 
!! 
!! @author J. Daniel Gezelter
!! @author Charles F. Vardeman II
!! @author Matthew Meineke
!! @version $Id: vector_class.F90,v 1.7 2005-04-15 22:03:59 gezelter Exp $, $Date: 2005-04-15 22:03:59 $, $Name: not supported by cvs2svn $, $Revision: 1.7 $

module Vector_class

  implicit NONE 
  PRIVATE

  public :: initialize
  public :: destroy
  public :: getSize
  public :: getElementAt
  public :: getPropertyListSize
  public :: getPropertyNameAt
  public :: addElement
  public :: setElementProperty
  public :: getElementProperty
  public :: getMatchingElementList
  public :: getFirstMatchingElement


  integer, parameter :: logical_data_type = 1
  integer, parameter :: integer_data_type = 2
  integer, parameter :: real_data_type = 3

  !! 
  type, public :: Vector 
     PRIVATE
     integer :: initialCapacity = 10
     integer :: capacityIncrement = 0
     integer :: elementCount = 0

     integer :: initialProperties = 5
     integer :: PropertyIncrement = 0
     integer :: propertyCount = 0

     integer, pointer :: ElementData(:) => null()
     character(len=100), pointer :: PropertyDescriptions(:) => null()
     integer, pointer :: PropertyDataType(:) => null()
     real(kind = 8), pointer :: realElementProperties(:,:) => null()
     integer, pointer :: integerElementProperties(:,:) => null()
     logical, pointer :: logicalElementProperties(:,:) => null()
  end type Vector

  !! Initialize vector
  interface initialize
     module procedure initialize_0i
     module procedure initialize_1i
     module procedure initialize_2i
     module procedure initialize_3i
     module procedure initialize_4i
  end interface

  interface setElementProperty
     module procedure setElementPropertyReal
     module procedure setElementPropertyInt
     module procedure setElementPropertyLogical
  end interface

  interface getElementProperty
     module procedure getElementPropertyReal
     module procedure getElementPropertyInt
     module procedure getElementPropertyLogical
  end interface

  interface getMatchingElementList
     module procedure getMatchingElementList1i
     module procedure getMatchingElementList2i
     module procedure getMatchingElementList1l
     module procedure getMatchingElementList2l
  end interface

  interface getFirstMatchingElement
     module procedure getFirstMatchingElement1i
     module procedure getFirstMatchingElement2i
     module procedure getFirstMatchingElement1l
     module procedure getFirstMatchingElement2l
  end interface
contains

  function getSize(this) result (ne)
    type(Vector), pointer :: this
    integer :: ne
    ne = this%elementCount
  end function getSize

  function getElementAt(this, loc) result (id)
    type(Vector), pointer :: this
    integer, intent(in) :: loc
    integer :: id
    id = this%ElementData(loc)
  end function getElementAt

  function getPropertyListSize(this) result (np)
    type(Vector), pointer :: this
    integer :: np
    np = this%propertyCount
  end function getPropertyListSize

  function getPropertyNameAt(this, loc) result (pn)
    type(Vector), pointer :: this
    integer, intent(in) :: loc
    character(len=len(this%PropertyDescriptions)) :: pn
    pn = this%PropertyDescriptions(loc)
  end function getPropertyNameAt

  function getFirstMatchingElement1i(this, MatchName, MatchValue) result (id)
    type(Vector), pointer :: this
    character(len=*), intent(in) :: MatchName
    integer, intent(in) :: MatchValue
    integer :: id
    integer :: i, j

    id = 0

    do i = 1, this%propertyCount
       if (this%PropertyDescriptions(i) == MatchName) then
          do j = 1, this%elementCount
             if (this%integerElementProperties(j, i) == MatchValue) then
                id = j
                return
             endif
          enddo
       endif
    enddo
    return
  end function getFirstMatchingElement1i

  function getFirstMatchingElement2i(this, MatchName1, MatchValue1, &
       MatchName2, MatchValue2) result (id)
    type(Vector), pointer :: this
    character(len=*), intent(in) :: MatchName1, MatchName2
    integer, intent(in) :: MatchValue1, MatchValue2
    integer :: id
    integer :: i, j, MatchID1, MatchID2
    logical :: found1 = .false.
    logical :: found2 = .false.

    id = 0
    ! first figure out which properties we are using to do the match:

    do i = 1, this%propertyCount
       if (this%PropertyDescriptions(i) == MatchName1) then
          MatchID1 = i
          found1 = .true.
       endif
       if (this%PropertyDescriptions(i) == MatchName2) then
          MatchID2 = i
          found2 = .true.
       endif

       if (found1.and.found2) then
          do j = 1, this%elementCount
             if ((this%integerElementProperties(j, MatchID1) == MatchValue1) &
                  .and. &
                  (this%integerElementProperties(j, MatchID2) ==MatchValue2)) &
                  then
                id = j
                return
             endif
          enddo
       endif
    end do

    return
  end function getFirstMatchingElement2i

  function getFirstMatchingElement1l(this, MatchName, MatchValue) result (id)
    type(Vector), pointer :: this
    character(len=*), intent(in) :: MatchName
    logical, intent(in) :: MatchValue
    integer :: id
    integer :: i, j

    id = 0

    do i = 1, this%propertyCount
       if (this%PropertyDescriptions(i) == MatchName) then
          do j = 1, this%elementCount
             if (this%logicalElementProperties(j, i) .eqv. MatchValue) then
                id = j
                return
             endif
          enddo
       endif
    enddo
    return
  end function getFirstMatchingElement1l

  function getFirstMatchingElement2l(this, MatchName1, MatchValue1, &
       MatchName2, MatchValue2) result (id)
    type(Vector), pointer :: this
    character(len=*), intent(in) :: MatchName1, MatchName2
    logical, intent(in) :: MatchValue1, MatchValue2
    integer :: id
    integer :: i, j, MatchID1, MatchID2
    logical :: found1 = .false.
    logical :: found2 = .false.

    id = 0
    ! first figure out which properties we are using to do the match:

    do i = 1, this%propertyCount
       if (this%PropertyDescriptions(i) == MatchName1) then
          MatchID1 = i
          found1 = .true.
       endif
       if (this%PropertyDescriptions(i) == MatchName2) then
          MatchID2 = i
          found2 = .true.
       endif

       if (found1.and.found2) then
          do j = 1, this%elementCount
             if ((this%logicalElementProperties(j, MatchID1).eqv.MatchValue1) &
                  .and. &
                  (this%logicalElementProperties(j, MatchID2).eqv.MatchValue2)) &
                  then
                id = j
                return
             endif
          enddo
       endif
    end do

    return
  end function getFirstMatchingElement2l

  subroutine getMatchingElementList1i(this, MatchName, MatchValue, &
       nMatches, MatchList)
    type(Vector), pointer :: this
    character(len=*), intent(in) :: MatchName
    integer, intent(in) :: MatchValue
    integer, intent(out) :: nMatches
    integer, pointer :: MatchList(:)
    integer :: i

    ! first figure out which property we are using to do the match:

    do i = 1, this%propertyCount
       if (this%PropertyDescriptions(i) == MatchName) then
          call getAllMatches1i(this, i, MatchValue, nMatches, MatchList)
          return
       endif
    enddo
    return
  end subroutine getMatchingElementList1i

  subroutine getMatchingElementList2i(this, MatchName1, MatchValue1, &
       MatchName2, MatchValue2, nMatches, MatchList)
    type(Vector), pointer :: this
    character(len=*), intent(in) :: MatchName1, MatchName2
    integer, intent(in)  :: MatchValue1, MatchValue2
    integer, intent(out)  :: nMatches
    integer, pointer :: MatchList(:)
    integer :: i, MatchID1, MatchID2
    logical :: found1 = .false.
    logical :: found2 = .false.

    ! first figure out which properties we are using to do the match:

    do i = 1, this%propertyCount
       if (this%PropertyDescriptions(i) == MatchName1) then
          MatchID1 = i
          found1 = .true.
       endif
       if (this%PropertyDescriptions(i) == MatchName2) then
          MatchID2 = i
          found2 = .true.
       endif

       if (found1.and.found2) then
          call getAllMatches2i(this, MatchID1, MatchValue1, &
               MatchID2, MatchValue2, nMatches, MatchList)
          return
       endif
    enddo
    return
  end subroutine getMatchingElementList2i

  subroutine getMatchingElementList1l(this, MatchName, MatchValue, &
       nMatches, MatchList)
    type(Vector), pointer :: this
    character(len=*), intent(in) :: MatchName
    logical, intent(in) :: MatchValue
    integer, intent(out) :: nMatches
    integer, pointer :: MatchList(:)
    integer :: i

    ! first figure out which property we are using to do the match:

    do i = 1, this%propertyCount
       if (this%PropertyDescriptions(i) == MatchName) then
          call getAllMatches1l(this, i, MatchValue, nMatches, MatchList)
          return
       endif
    enddo
    return
  end subroutine getMatchingElementList1l

  subroutine getMatchingElementList2l(this, MatchName1, MatchValue1, &
       MatchName2, MatchValue2, nMatches, MatchList)
    type(Vector), pointer :: this
    character(len=*), intent(in) :: MatchName1, MatchName2
    logical, intent(in)  :: MatchValue1, MatchValue2
    integer, intent(out)  :: nMatches
    integer, pointer :: MatchList(:)
    integer :: i, MatchID1, MatchID2
    logical :: found1 = .false.
    logical :: found2 = .false.

    ! first figure out which properties we are using to do the match:

    do i = 1, this%propertyCount
       if (this%PropertyDescriptions(i) == MatchName1) then
          MatchID1 = i
          found1 = .true.
       endif
       if (this%PropertyDescriptions(i) == MatchName2) then
          MatchID2 = i
          found2 = .true.
       endif

       if (found1.and.found2) then
          call getAllMatches2l(this, MatchID1, MatchValue1, &
               MatchID2, MatchValue2, nMatches, MatchList)
          return
       endif
    enddo
    return
  end subroutine getMatchingElementList2l

  subroutine getAllMatches1i(this, MatchID, MatchValue, nMatches, MatchList)
    type(Vector), pointer :: this
    integer, intent(in) :: MatchID
    integer, intent(in) :: MatchValue
    integer, pointer :: MatchList(:) 
    integer, allocatable :: MatchListTemp(:)
    integer, intent(out) :: nMatches
    integer :: error, i

    if(associated(MatchList)) deallocate(MatchList)
    MatchList => null()

    allocate(MatchListTemp(this%elementCount), stat=error)
    if(error .ne. 0) write(*,*) 'Could not allocate MatchListTemp!'

    nMatches = 0

    do i = 1, this%elementCount
       if (this%integerElementProperties(i, MatchID) == MatchValue) then
          nMatches = nMatches + 1
          MatchListTemp(nMatches) = i
       endif
    enddo


    if (nMatches .ne. 0) then
       allocate(MatchList(nMatches), stat=error)
       if (error.ne.0) write(*, *) 'Could not allocate MatchList!'
       do i = 1, nMatches
          MatchList(i) = MatchListTemp(i)
       enddo
    endif

    deallocate(MatchListTemp)


  end subroutine getAllMatches1i

  subroutine getAllMatches2i(this, MatchID1, MatchValue1, &
       MatchID2, MatchValue2, nMatches, MatchList)
    type(Vector), pointer :: this
    integer, intent(in) :: MatchID1, MatchID2
    integer, intent(in) :: MatchValue1, MatchValue2
    integer, pointer :: MatchList(:)
    integer, allocatable :: MatchListTemp(:)
    integer, intent(out) :: nMatches
    integer :: error, i

    if(associated(MatchList)) deallocate(MatchList)
    MatchList => null()

    allocate(MatchListTemp(this%elementCount), stat=error)
    if(error .ne. 0) write(*,*) 'Could not allocate MatchListTemp!'

    nMatches = 0

    do i = 1, this%elementCount
       if ((this%integerElementProperties(i, MatchID1) == MatchValue1) .and. &
            (this%integerElementProperties(i, MatchID2) == MatchValue2)) then
          nMatches = nMatches + 1
          MatchListTemp(nMatches) = i
       endif
    enddo

    if (nMatches .ne. 0) then
       allocate(MatchList(nMatches), stat=error)
       if (error.ne.0) write(*, *) 'Could not allocate MatchList!'
       do i = 1, nMatches
          MatchList(i) = MatchListTemp(i)
       enddo
    endif

    deallocate(MatchListTemp)

  end subroutine getAllMatches2i

  subroutine getAllMatches1l(this, MatchID, MatchValue, nMatches, MatchList)
    type(Vector), pointer :: this
    integer, intent(in) :: MatchID
    logical, intent(in) :: MatchValue
    integer, pointer :: MatchList(:)
    integer, allocatable :: MatchListTemp(:)
    integer, intent(out) :: nMatches
    integer :: error, i

    if(associated(MatchList)) deallocate(MatchList)
    MatchList => null()

    allocate(MatchListTemp(this%elementCount), stat=error)
    if(error .ne. 0) write(*,*) 'Could not allocate MatchListTemp!'

    nMatches = 0

    do i = 1, this%elementCount
       if (this%logicalElementProperties(i, MatchID).eqv.MatchValue) then
          nMatches = nMatches + 1
          MatchListTemp(nMatches) = i
       endif
    enddo

    if (nMatches .ne. 0) then
       allocate(MatchList(nMatches), stat=error)
       if (error.ne.0) write(*, *) 'Could not allocate MatchList!'
       do i = 1, nMatches
          MatchList(i) = MatchListTemp(i)
       enddo
    endif

    deallocate(MatchListTemp)

  end subroutine getAllMatches1l

  subroutine getAllMatches2l(this, MatchID1, MatchValue1, &
       MatchID2, MatchValue2, nMatches, MatchList)
    type(Vector), pointer :: this
    integer, intent(in) :: MatchID1, MatchID2
    logical, intent(in) :: MatchValue1, MatchValue2
    integer, pointer :: MatchList(:)
    integer, allocatable :: MatchListTemp(:)
    integer, intent(out) :: nMatches
    integer :: error, i

    if(associated(MatchList)) deallocate(MatchList)
    MatchList => null()

    allocate(MatchListTemp(this%elementCount), stat=error)
    if(error .ne. 0) write(*,*) 'Could not allocate MatchListTemp!'

    nMatches = 0

    do i = 1, this%elementCount
       if ((this%logicalElementProperties(i, MatchID1).eqv.MatchValue1) .and. &
            (this%logicalElementProperties(i, MatchID2).eqv.MatchValue2)) then
          nMatches = nMatches + 1
          MatchListTemp(nMatches) = i
       endif
    enddo

    if (nMatches .ne. 0) then
       allocate(MatchList(nMatches), stat=error)
       if (error.ne.0) write(*, *) 'Could not allocate MatchList!'
       do i = 1, nMatches
          MatchList(i) = MatchListTemp(i)
       enddo
    endif

    deallocate(MatchListTemp)

  end subroutine getAllMatches2l


  subroutine getElementPropertyReal(this, id, PropName, pv)
    type(Vector), pointer :: this
    integer :: id, whichprop
    character(len=*) :: PropName
    real( kind = 8 ) :: pv

    whichprop = getPropertyIndex(this, PropName)
    if (whichprop .eq. 0 ) then 
       write(*,*) 'unknown property: ', PropName
       pv = 0.0
    else
       if (this%PropertyDataType(whichprop) .ne. real_data_type) then
          write(*,*) 'Property: ', PropName, " is not real data type."
          pv = 0.0
       else
          pv = this%realElementProperties(id, whichprop)
       endif
    endif
  end subroutine getElementPropertyReal

  subroutine getElementPropertyInt(this, id, PropName, pv)
    type(Vector), pointer :: this
    integer :: id, whichprop
    character(len=*) :: PropName
    integer :: pv

    whichprop = getPropertyIndex(this, PropName)
    if (whichprop .eq. 0 ) then 
       write(*,*) 'unknown property! ', PropName
       pv = 0
    else
       if (this%PropertyDataType(whichprop) .ne. integer_data_type) then
          write(*,*) 'Property! ', PropName, " is not integer data type."
          pv = 0
       else
          pv = this%integerElementProperties(id, whichprop)
       endif
    endif
  end subroutine getElementPropertyInt

  subroutine getElementPropertyLogical(this, id, PropName, pv)
    type(Vector), pointer :: this
    integer :: id, whichprop
    character(len=*) :: PropName
    logical :: pv

    whichprop = getPropertyIndex(this, PropName)
    if (whichprop .eq. 0 ) then 
       write(*,*) 'unknown property! ', PropName
       pv = .false.
    else
       if (this%PropertyDataType(whichprop) .ne. logical_data_type) then
          write(*,*) 'Property! ', PropName, " is not logical data type."
          pv = .false.
       else
          pv = this%logicalElementProperties(id, whichprop)
       endif
    endif
  end subroutine getElementPropertyLogical

  function getPropertyIndex(this, PropName) result (id)
    type(Vector), pointer :: this
    integer :: id, i
    character(len=*) :: PropName

    do i = 1, this%propertyCount
       if (this%PropertyDescriptions(i) == PropName) then
          id = i
          return
       endif
    enddo
    id = 0
  end function getPropertyIndex

  subroutine ensureCapacityHelper(this, minCapacity, minPropCap)
    type(Vector), pointer :: this, that
    integer, intent(in) :: minCapacity, minPropCap
    integer :: oldCapacity, oldPropCap 
    integer :: newCapacity, newPropCap
    logical :: resizeFlag 

    resizeFlag = .false.

    !  first time: allocate a new vector with default size

    if (.not. associated(this)) then
       this => initialize()
    endif

    oldCapacity = size(this%ElementData)
    oldPropCap  = size(this%PropertyDescriptions)

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

!!! newCapacity is not set.....
    if (minPropCap > oldPropCap) then
       if (this%PropertyIncrement .gt. 0) then
          newPropCap = oldPropCap + this%PropertyIncrement
       else
          newPropCap = oldPropCap * 2
       endif
       if (newPropCap .lt. minPropCap) then
          newPropCap = minPropCap
       endif
       resizeFlag = .true.
    else
       newPropCap = oldPropCap
    endif

    if (resizeFlag) then
       that => initialize(newCapacity, newPropCap, &
            this%capacityIncrement, this%PropertyIncrement)       
       call copyAllData(this, that)
       this => destroy(this)
       this => that
    endif
  end subroutine ensureCapacityHelper

  subroutine copyAllData(v1, v2)
    type(Vector), pointer  :: v1
    type(Vector), pointer  :: v2
    integer :: i, j

    do i = 1, v1%elementCount
       v2%elementData(i) = v1%elementData(i)
       do j = 1, v1%propertyCount

          if (v1%PropertyDataType(j) .eq. integer_data_type) &
               v2%integerElementProperties(i,j) = &
               v1%integerElementProperties(i,j)

          if (v1%PropertyDataType(j) .eq. real_data_type)  &
               v2%realElementProperties(i,j) = v1%realElementProperties(i,j)

          if (v1%PropertyDataType(j) .eq. logical_data_type)  &
               v2%logicalElementProperties(i,j) = &
               v1%logicalElementProperties(i,j)                        
       enddo
    enddo

    do j = 1, v1%propertyCount
       v2%PropertyDescriptions(j) = v1%PropertyDescriptions(j)
       v2%PropertyDataType(j) = v1%PropertyDataType(j)
    enddo

    v2%elementCount = v1%elementCount
    v2%propertyCount = v1%propertyCount

    return
  end subroutine copyAllData

  function addElement(this) result (id)
    type(Vector), pointer :: this
    integer :: id
    integer :: error

    if (.not. associated(this)) then
       call ensureCapacityHelper(this,1,0)
    else
       call ensureCapacityHelper(this, this%elementCount + 1, this%PropertyCount)
    end if

    this%elementCount = this%elementCount + 1

    !! We never use this and we set the entire array to the same value
    this%elementData = this%elementCount
    id = this%elementCount
  end function addElement

  recursive subroutine setElementPropertyReal(this, id, PropName, PropValue)
    type(Vector), pointer :: this
    integer :: id, i
    character(len=*), intent(in) :: PropName
    real( kind = 8 ), intent(in) :: PropValue
    logical :: foundit

    foundit = .false.

    ! first make sure that the PropName isn't in the list of known properties: 

    do i = 1, this%propertyCount
       if (PropName == this%PropertyDescriptions(i)) then
          foundit = .true.
          this%realElementProperties(id,i) = PropValue
       endif
    enddo

    if (.not.foundit) then
       call addPropertyToVector(this, PropName, real_data_type)
       call setElementPropertyReal(this, id, PropName, PropValue)
    endif
  end subroutine setElementPropertyReal

  recursive subroutine setElementPropertyInt(this, id, PropName, PropValue)
    type(Vector), pointer :: this
    integer :: id, i
    character(len=*), intent(in) :: PropName
    integer, intent(in) :: PropValue
    logical :: foundit

    foundit = .false.
    ! first make sure that the PropName isn't in the list of known properties: 
    do i = 1, this%propertyCount
       if (PropName == this%PropertyDescriptions(i)) then
          foundit = .true.
          this%integerElementProperties(id,i) = PropValue
       endif
    enddo

    if (.not.foundit) then
       call addPropertyToVector(this, PropName, integer_data_type)
       call setElementPropertyInt(this, id, PropName, PropValue)
    endif
  end subroutine setElementPropertyInt

  recursive subroutine setElementPropertyLogical(this, id, PropName, PropValue)
    type(Vector), pointer :: this
    integer :: id, i
    character(len=*), intent(in) :: PropName
    logical, intent(in) :: PropValue
    logical :: foundit 

    foundit = .false.
    ! first make sure that the PropName isn't in the list of known properties: 
    do i = 1, this%propertyCount
       if (PropName == this%PropertyDescriptions(i)) then
          foundit = .true.
          this%logicalElementProperties(id,i) = PropValue
       endif
    enddo

    if (.not.foundit) then
       call addPropertyToVector(this, PropName, logical_data_type)
       call setElementPropertyLogical(this, id, PropName, PropValue)
    endif
  end subroutine setElementPropertyLogical

  subroutine addPropertyToVector(this, PropName, data_type)
    type(Vector), pointer :: this
    character(len=*), intent(in) :: PropName
    integer data_type

    call ensureCapacityHelper(this, this%elementCount, this%propertyCount + 1)
    this%propertyCount = this%propertyCount + 1
    this%PropertyDescriptions(this%propertyCount) = PropName
    this%PropertyDataType(this%propertyCount) = data_type
  end subroutine addPropertyToVector

  function initialize_0i() result(this)
    type(Vector), pointer :: this 
    this => initialize_2i(10, 5)
  end function initialize_0i

  function initialize_1i(nprop) result(this)
    integer, intent(in) :: nprop
    type(Vector), pointer :: this 
    this => initialize_2i(10, nprop)
  end function initialize_1i

  function initialize_2i(cap, nprop) result(this)
    integer, intent(in) :: cap, nprop
    type(Vector), pointer :: this 
    this => initialize_4i(cap, nprop, 0, 0)
  end function initialize_2i

  function initialize_3i(cap, nprop, capinc) result(this)
    integer, intent(in) :: cap, nprop, capinc
    type(Vector), pointer :: this 
    this => initialize_4i(cap, nprop, capinc, 0)
  end function initialize_3i

  function initialize_4i(cap, nprop, capinc, propinc) result(this)
    integer, intent(in) :: cap, nprop, capinc, propinc
    integer :: error
    type(Vector), pointer :: this 

    nullify(this)

    if (cap < 0) then
       write(*,*) 'Bogus Capacity:', cap
       return
    endif
    if (nprop < 0) then
       write(*,*) 'Bogus Number of Properties:', nprop
       return
    endif

    allocate(this,stat=error)
    if ( error /= 0 ) then
       write(*,*) 'Could not allocate Vector!'
       return
    end if

    this%initialCapacity = cap
    this%initialProperties = nprop
    this%capacityIncrement = capinc
    this%propertyIncrement = propinc

    allocate(this%elementData(this%initialCapacity), stat=error)
    if(error /= 0) write(*,*) 'Could not allocate elementData!'


    allocate(this%PropertyDescriptions(this%initialProperties), &
         stat=error)
    if(error /=  0) write(*,*) 'Could not allocate PropertyDescriptions!'

    allocate(this%PropertyDataType(this%initialProperties), &
         stat=error)
    if(error /=  0) write(*,*) 'Could not allocate PropertyDataType!'

    allocate(this%integerElementProperties(this%initialCapacity, &
         this%initialProperties), stat=error)
    if(error /= 0) write(*,*) 'Could not allocate integerElementProperties!'

    allocate(this%realElementProperties(this%initialCapacity, &
         this%initialProperties), stat=error)
    if(error /= 0) write(*,*) 'Could not allocate realElementProperties!'   

    allocate(this%logicalElementProperties(this%initialCapacity, &
         this%initialProperties), stat=error)
    if(error .ne. 0) write(*,*) 'Could not allocate logicalElementProperties!'

  end function initialize_4i

  !! This function destroys the vector components....
  function destroy(this) result(null_this)
    logical :: done
    type(Vector), pointer :: this 
    type(Vector), pointer :: null_this 

    if (.not. associated(this)) then
       null_this => null()
       return
    end if

    !! Walk down the list and deallocate each of the vector component
    if(associated(this%logicalElementProperties)) then
       deallocate(this%logicalElementProperties)
       this%logicalElementProperties=>null()
    endif
    if(associated(this%realElementProperties)) then
       deallocate(this%realElementProperties)
       this%realElementProperties=>null()
    endif
    if(associated(this%integerElementProperties)) then
       deallocate(this%integerElementProperties)
       this%integerElementProperties=>null()
    endif
    if(associated(this%PropertyDataType)) then
       deallocate(this%PropertyDataType)
       this%PropertyDataType=>null()
    endif
    if(associated(this%PropertyDescriptions)) then
       deallocate(this%PropertyDescriptions)
       this%PropertyDescriptions=>null()
    endif
    if(associated(this%elementData)) then
       deallocate(this%elementData)
       this%elementData=>null()
    endif
    deallocate(this)
    this => null()
    null_this => null()
  end function destroy

end module Vector_class
