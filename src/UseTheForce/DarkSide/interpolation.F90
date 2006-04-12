!!
!! Copyright (c) 2006 The University of Notre Dame. All Rights Reserved.
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
!!
!!  interpolation.F90
!!
!!  Created by Charles F. Vardeman II on 03 Apr 2006.
!!
!!  PURPOSE: Generic Spline interplelation routines. These routines assume that we are on a uniform grid for
!!           precomputation of spline parameters.
!!
!! @author Charles F. Vardeman II 
!! @version $Id: interpolation.F90,v 1.1 2006-04-12 21:15:17 chuckv Exp $


module  INTERPOLATION
  use definitions
  use status
  implicit none
  PRIVATE

  character(len = statusMsgSize) :: errMSG

  type, public :: splineType
     private
     integer :: npoints = 0
     real(kind=dp) :: delta_x
     real(kind=dp) :: range
     real(kind=dp) :: range_inv
     real (kind=dp), pointer,dimension(:) :: xa => null()
     real (kind=dp), pointer,dimension(:) :: ya => null()
     real (kind=dp), pointer,dimension(:) :: yppa => null()
  end type splineType

  type, public :: multiSplineType
     private
     integer :: npoints = 0
     integer :: nfuncs = 0

     integer :: npoints = 0
     real(kind=dp) :: delta_x
     real(kind=dp) :: range
     real(kind=dp) :: range_inv
     real (kind=dp), pointer,dimension(:)   :: xa => null()
     real (kind=dp), pointer,dimension(:,:) :: ya => null()
     real (kind=dp), pointer,dimension(:,:) :: yppa => null()
  end type splineType


  interface splineLookup
     module procedure multiSplint
     module procedure splintd
     module procedure splintd1
     module procedure splintd2
  end interface

  public :: splint
  public :: newSpline
  public :: newMultiSpline
  public :: deleteSpline
  public :: deleteMultiSpline
  

contains

  !! mySpline is pointer to spline type, nx number of data points, 
  !! xa tabulated x values and ya respective values for xa, yp1 
  !! is value for derivative at first point and ypn is value 
  !! for derivative at point n.
  subroutine newSpline(thisSpline,nx, xa, ya, yp1, ypn, boundary)

    ! yp1 and ypn are the first derivatives of y at the two endpoints
    ! if boundary is 'L' the lower derivative is used
    ! if boundary is 'U' the upper derivative is used
    ! if boundary is 'B' then both derivatives are used
    ! if boundary is anything else, then both derivatives are assumed to be 0


    type (splineType), intent(inout) :: thisSpline


    real( kind = DP ), pointer, dimension(:)        :: xa
    real( kind = DP ), pointer, dimension(:)        :: ya
    real( kind = DP ), dimension(size(xa)) :: u
    real( kind = DP ) :: yp1,ypn,un,qn,sig,p
    character(len=*) :: boundary
    integer :: nx, i, k, max_array_size
    integer :: alloc_error

    alloc_error = 0

    if (thisSpline%npoints/=0) then
       call handleWarning("INTERPOLATION:newSpline",&
            "Type has already been created")
       call deleteSpline(thisSpline)
    end if


    ! make sure the sizes match
    if ((nx /= size(xa)) .or. (nx /= size(ya))) then
       call handleWarning("INTERPOLATION:newSpline","Array size mismatch")
    end if

   
    thisSpline%npoints = nx
    allocate(thisSpline%yppa(nx),stat=alloc_error)
    if(alloc_error .ne. 0) call handleWarning("INTERPOLATION:newSpline",&
         "Error in allocating storage for yppa")

    thisSpline%xa => xa
    thisSpline%ya => ya




    if ((boundary.eq.'l').or.(boundary.eq.'L').or. &
         (boundary.eq.'b').or.(boundary.eq.'B')) then
       thisSpline%yppa(1) = -0.5E0_DP
       u(1) = (3.0E0_DP/(xa(2)-xa(1)))*((ya(2)-&
            ya(1))/(xa(2)-xa(1))-yp1)
    else
       thisSpline%yppa(1) = 0.0E0_DP
       u(1)  = 0.0E0_DP
    endif

    do i = 2, nx - 1
       sig = (thisSpline%xa(i) - thisSpline%xa(i-1)) / (thisSpline%xa(i+1) - thisSpline%xa(i-1))
       p = sig * thisSpline%yppa(i-1) + 2.0E0_DP
       thisSpline%yppa(i) = (sig - 1.0E0_DP) / p
       u(i) = (6.0E0_DP*((thisSpline%ya(i+1)-thisSpline%ya(i))/(thisSpline%xa(i+1)-thisSpline%xa(i)) - &
            (thisSpline%ya(i)-thisSpline%ya(i-1))/(thisSpline%xa(i)-thisSpline%xa(i-1)))/ &
            (thisSpline%xa(i+1)-thisSpline%xa(i-1)) - sig * u(i-1))/p
    enddo

    if ((boundary.eq.'u').or.(boundary.eq.'U').or. &
         (boundary.eq.'b').or.(boundary.eq.'B')) then
       qn = 0.5E0_DP
       un = (3.0E0_DP/(thisSpline%xa(nx)-thisSpline%xa(nx-1)))* &
            (ypn-(thisSpline%ya(nx)-thisSpline%ya(nx-1))/(thisSpline%xa(nx)-thisSpline%xa(nx-1)))
    else
       qn = 0.0E0_DP
       un = 0.0E0_DP
    endif

    thisSpline%yppa(nx)=(un-qn*u(nx-1))/(qn*thisSpline%yppa(nx-1)+1.0E0_DP)

    do k = nx-1, 1, -1
       thisSpline%yppa(k)=thisSpline%yppa(k)*thisSpline%yppa(k+1)+u(k)
    enddo

  end subroutine newSpline

  subroutine deleteSpline(thisSpline)
    type(splineType) :: thisSpline

  
    
    if(associated(thisSpline%xa)) then
       deallocate(thisSpline%xa)
       thisSpline%xa => null()
    end if
    if(associated(thisSpline%ya)) then
       deallocate(thisSpline%ya)
       thisSpline%ya => null()
    end if
    if(associated(thisSpline%yppa)) then
       deallocate(thisSpline%yppa)
       thisSpline%yppa => null()
    end if    
    
    thisSpline%npoints=0
    
  end subroutine deleteSpline

   subroutine splintd2(thisSpline, x, y, dy, d2y)
    type(splineType) :: thisSpline
    real( kind = DP ), intent(in) :: x
    real( kind = DP ), intent(out) :: y,dy,d2y

    
    real( kind = DP ) :: del, h, a, b, c, d
    integer :: j

    ! this spline code assumes that the x points are equally spaced
    ! do not attempt to use this code if they are not.


    ! find the closest point with a value below our own:
    j = FLOOR(real((thisSpline%npoints-1),kind=dp) * &
         (x - thisSpline%xa(1)) / (thisSpline%xa(thisSpline%npoints) - thisSpline%xa(1))) + 1

    ! check to make sure we're inside the spline range:
    if ((j.gt.thisSpline%npoints).or.(j.lt.1)) then
       write(errMSG,*) "EAM_splint: x is outside bounds of spline: ",x,j
       call handleError("INTERPOLATION::SPLINT2d",errMSG)
    endif
    ! check to make sure we haven't screwed up the calculation of j:
    if ((x.lt.thisSpline%xa(j)).or.(x.gt.thisSpline%xa(j+1))) then
       if (j.ne.thisSpline%npoints) then
          write(errMSG,*) "EAM_splint:",x," x is outside bounding range"
          call handleError("INTERPOLATION::SPLINT2d",errMSG)
       endif
    endif

    del = thisSpline%xa(j+1) - x
    h = thisSpline%xa(j+1) - thisSpline%xa(j)

    a = del / h
    b = 1.0E0_DP - a
    c = a*(a*a - 1.0E0_DP)*h*h/6.0E0_DP
    d = b*(b*b - 1.0E0_DP)*h*h/6.0E0_DP

    y = a*thisSpline%ya(j) + b*thisSpline%ya(j+1) + c*thisSpline%yppa(j) + d*thisSpline%yppa(j+1)

    dy = (thisSpline%ya(j+1)-thisSpline%ya(j))/h &
         - (3.0E0_DP*a*a - 1.0E0_DP)*h*thisSpline%yppa(j)/6.0E0_DP &
         + (3.0E0_DP*b*b - 1.0E0_DP)*h*thisSpline%yppa(j+1)/6.0E0_DP


    d2y = a*thisSpline%yppa(j) + b*thisSpline%yppa(j+1)


  end subroutine splintd2
   subroutine splintd1(thisSpline, x, y, dy)
    type(splineType) :: thisSpline
    real( kind = DP ), intent(in) :: x
    real( kind = DP ), intent(out) :: y,dy

    
    real( kind = DP ) :: del, h, a, b, c, d
    integer :: j

    ! this spline code assumes that the x points are equally spaced
    ! do not attempt to use this code if they are not.


    ! find the closest point with a value below our own:
    j = FLOOR(real((thisSpline%npoints-1),kind=dp) *&
         (x - thisSpline%xa(1)) / (thisSpline%xa(thisSpline%npoints) - thisSpline%xa(1))) + 1

    ! check to make sure we're inside the spline range:
    if ((j.gt.thisSpline%npoints).or.(j.lt.1)) then
       write(errMSG,*) "EAM_splint: x is outside bounds of spline: ",x,j
       call handleError("INTERPOLATION::SPLINT2d",errMSG)
    endif
    ! check to make sure we haven't screwed up the calculation of j:
    if ((x.lt.thisSpline%xa(j)).or.(x.gt.thisSpline%xa(j+1))) then
       if (j.ne.thisSpline%npoints) then
          write(errMSG,*) "EAM_splint:",x," x is outside bounding range"
          call handleError("INTERPOLATION::SPLINT2d",errMSG)
       endif
    endif

    del = thisSpline%xa(j+1) - x
    h = thisSpline%xa(j+1) - thisSpline%xa(j)

    a = del / h
    b = 1.0E0_DP - a
    c = a*(a*a - 1.0E0_DP)*h*h/6.0E0_DP
    d = b*(b*b - 1.0E0_DP)*h*h/6.0E0_DP

    y = a*thisSpline%ya(j) + b*thisSpline%ya(j+1) + c*thisSpline%yppa(j) + d*thisSpline%yppa(j+1)

    dy = (thisSpline%ya(j+1)-thisSpline%ya(j))/h &
         - (3.0E0_DP*a*a - 1.0E0_DP)*h*thisSpline%yppa(j)/6.0E0_DP &
         + (3.0E0_DP*b*b - 1.0E0_DP)*h*thisSpline%yppa(j+1)/6.0E0_DP


  


  end subroutine splintd1
   subroutine splintd(thisSpline, x, y)
    type(splineType) :: thisSpline
    real( kind = DP ), intent(in) :: x
    real( kind = DP ), intent(out) :: y

    
    real( kind = DP ) :: del, h, a, b, c, d
    integer :: j

    ! this spline code assumes that the x points are equally spaced
    ! do not attempt to use this code if they are not.


    ! find the closest point with a value below our own:
    j = FLOOR(real((thisSpline%npoints-1),kind=dp) * &
         (x - thisSpline%xa(1)) / (thisSpline%xa(thisSpline%npoints) - thisSpline%xa(1))) + 1

    ! check to make sure we're inside the spline range:
    if ((j.gt.thisSpline%npoints).or.(j.lt.1)) then
       write(errMSG,*) "EAM_splint: x is outside bounds of spline: ",x,j
       call handleError("INTERPOLATION::SPLINT2d",errMSG)
    endif
    ! check to make sure we haven't screwed up the calculation of j:
    if ((x.lt.thisSpline%xa(j)).or.(x.gt.thisSpline%xa(j+1))) then
       if (j.ne.thisSpline%npoints) then
          write(errMSG,*) "EAM_splint:",x," x is outside bounding range"
          call handleError("INTERPOLATION::SPLINT2d",errMSG)
       endif
    endif

    del = thisSpline%xa(j+1) - x
    h = thisSpline%xa(j+1) - thisSpline%xa(j)

    a = del / h
    b = 1.0E0_DP - a
    c = a*(a*a - 1.0E0_DP)*h*h/6.0E0_DP
    d = b*(b*b - 1.0E0_DP)*h*h/6.0E0_DP

    y = a*thisSpline%ya(j) + b*thisSpline%ya(j+1) + c*thisSpline%yppa(j) + d*thisSpline%yppa(j+1)

  end subroutine splintd
  

end module INTERPOLATION
