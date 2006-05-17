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
!!  PURPOSE: Generic Spline interpolation routines. 
!!
!! @author Charles F. Vardeman II 
!! @version $Id: interpolation.F90,v 1.8 2006-05-17 15:37:15 gezelter Exp $


module interpolation
  use definitions
  use status
  implicit none
  PRIVATE

  type, public :: cubicSpline
     logical :: isUniform = .false.
     integer :: n = 0
     real(kind=dp) :: dx_i
     real (kind=dp), pointer,dimension(:)   :: x => null()
     real (kind=dp), pointer,dimension(:)   :: y => null()
     real (kind=dp), pointer,dimension(:)   :: b => null()
     real (kind=dp), pointer,dimension(:)   :: c => null()
     real (kind=dp), pointer,dimension(:)   :: d => null()
  end type cubicSpline

  public :: newSpline
  public :: deleteSpline
  public :: lookupSpline
  public :: lookupUniformSpline
  public :: lookupNonuniformSpline
  public :: lookupUniformSpline1d
  
contains
  

  subroutine newSpline(cs, x, y, isUniform)
    
    implicit none

    type (cubicSpline), intent(inout) :: cs
    real( kind = DP ), intent(in) :: x(:), y(:)
    real( kind = DP ) :: fp1, fpn, p
    REAL( KIND = DP), DIMENSION(size(x)-1) :: diff_y, H

    logical, intent(in) :: isUniform
    integer :: i, alloc_error, n, k

    alloc_error = 0

    if (cs%n .ne. 0) then
       call handleWarning("interpolation::newSpline", &
            "cubicSpline struct was already created")
       call deleteSpline(cs)
    end if

    ! make sure the sizes match

    n = size(x)
    
    if ( size(y) .ne. size(x) ) then
       call handleError("interpolation::newSpline", &
            "Array size mismatch")
    end if
    
    cs%n = n
    cs%isUniform = isUniform

    allocate(cs%x(n), stat=alloc_error)
    if(alloc_error .ne. 0) then
       call handleError("interpolation::newSpline", &
            "Error in allocating storage for x")
    endif

    allocate(cs%y(n), stat=alloc_error)
    if(alloc_error .ne. 0) then
       call handleError("interpolation::newSpline", &
            "Error in allocating storage for y")
    endif

    allocate(cs%b(n), stat=alloc_error)
    if(alloc_error .ne. 0) then
       call handleError("interpolation::newSpline", &
            "Error in allocating storage for b")
    endif

    allocate(cs%c(n), stat=alloc_error)
    if(alloc_error .ne. 0) then
       call handleError("interpolation::newSpline", &
            "Error in allocating storage for c")
    endif

    allocate(cs%d(n), stat=alloc_error)
    if(alloc_error .ne. 0) then
       call handleError("interpolation::newSpline", &
            "Error in allocating storage for d")
    endif

    ! make sure we are monotinically increasing in x:

    h = diff(x)
    if (any(h <= 0)) then
       call handleError("interpolation::newSpline", &
            "Negative dx interval found")
    end if

    ! load x and y values into the cubicSpline structure:

    do i = 1, n
       cs%x(i) = x(i)
       cs%y(i) = y(i)
    end do

    ! Calculate coefficients for the tridiagonal system: store
    ! sub-diagonal in B, diagonal in D, difference quotient in C.

    cs%b(1:n-1) = h
    diff_y = diff(y)
    cs%c(1:n-1) = diff_y / h

    if (n == 2) then
       ! Assume the derivatives at both endpoints are zero
       ! another assumption could be made to have a linear interpolant
       ! between the two points.  In that case, the b coefficients
       ! below would be diff_y(1)/h(1) and the c and d coefficients would
       ! both be zero.
       cs%b(1) = 0.0_dp
       cs%c(1) = -3.0_dp * (diff_y(1)/h(1))**2
       cs%d(1) = -2.0_dp * (diff_y(1)/h(1))**3
       cs%b(2) = cs%b(1)
       cs%c(2) = 0.0_dp
       cs%d(2) = 0.0_dp
       cs%dx_i = 1.0_dp / h(1)
      return
    end if

    cs%d(1) = 2.0_dp * cs%b(1)
    do i = 2, n-1
      cs%d(i) = 2.0_dp * (cs%b(i) + cs%b(i-1))
    end do
    cs%d(n) = 2.0_dp * cs%b(n-1)

    ! Calculate estimates for the end slopes using polynomials
    ! that interpolate the data nearest the end.
    
    fp1 = cs%c(1) - cs%b(1)*(cs%c(2) - cs%c(1))/(cs%b(1) + cs%b(2))
    if (n > 3) then
      fp1 = fp1 + cs%b(1)*((cs%b(1) + cs%b(2))*(cs%c(3) - cs%c(2))/ &
           (cs%b(2) + cs%b(3)) - cs%c(2) + cs%c(1))/(x(4) - x(1))
    end if
          
    fpn = cs%c(n-1) + cs%b(n-1)*(cs%c(n-1) - cs%c(n-2))/(cs%b(n-2) + cs%b(n-1))
    if (n > 3) then
      fpn = fpn + cs%b(n-1)*(cs%c(n-1) - cs%c(n-2) - (cs%b(n-2) + cs%b(n-1))* &
           (cs%c(n-2) - cs%c(n-3))/(cs%b(n-2) + cs%b(n-3)))/(x(n) - x(n-3))
    end if

    ! Calculate the right hand side and store it in C.

    cs%c(n) = 3.0_dp * (fpn - cs%c(n-1))
    do i = n-1,2,-1
      cs%c(i) = 3.0_dp * (cs%c(i) - cs%c(i-1))
    end do
    cs%c(1) = 3.0_dp * (cs%c(1) - fp1)

    ! Solve the tridiagonal system.

    do k = 2, n
      p = cs%b(k-1) / cs%d(k-1)
      cs%d(k) = cs%d(k) - p*cs%b(k-1)
      cs%c(k) = cs%c(k) - p*cs%c(k-1)
    end do
    cs%c(n) = cs%c(n) / cs%d(n)
    do k = n-1, 1, -1
      cs%c(k) = (cs%c(k) - cs%b(k) * cs%c(k+1)) / cs%d(k)
    end do

    ! Calculate the coefficients defining the spline.

    cs%d(1:n-1) = diff(cs%c) / (3.0_dp * h)
    cs%b(1:n-1) = diff_y / h - h * (cs%c(1:n-1) + h * cs%d(1:n-1))
    cs%b(n) = cs%b(n-1) + h(n-1) * (2.0_dp*cs%c(n-1) + h(n-1)*3.0_dp*cs%d(n-1))

    if (isUniform) then
       cs%dx_i = 1.0_dp / (x(2) - x(1))
    endif

    return
    
  contains
    
    function diff(v)
      ! Auxiliary function to compute the forward difference
      ! of data stored in a vector v.
      
      implicit none
      real (kind = dp), dimension(:), intent(in) :: v
      real (kind = dp), dimension(size(v)-1) :: diff
      
      integer :: n
      
      n = size(v)
      diff = v(2:n) - v(1:n-1)
      return
    end function diff
    
  end subroutine newSpline
       
  subroutine deleteSpline(this)

    type(cubicSpline) :: this
    
    if(associated(this%x)) then
       deallocate(this%x)
       this%x => null()
    end if
    if(associated(this%c)) then
       deallocate(this%c)
       this%c => null()
    end if
    
    this%n = 0
    
  end subroutine deleteSpline

  subroutine lookupNonuniformSpline(cs, xval, yval)
    
    implicit none

    type (cubicSpline), intent(in) :: cs
    real( kind = DP ), intent(in)  :: xval
    real( kind = DP ), intent(out) :: yval
    real( kind = DP ) ::  dx
    integer :: i, j
    !
    !  Find the interval J = [ cs%x(J), cs%x(J+1) ] that contains 
    !  or is nearest to xval.
    !
    j = cs%n - 1

    do i = 0, cs%n - 2

       if ( xval < cs%x(i+1) ) then
          j = i
          exit
       end if

    end do
    !
    !  Evaluate the cubic polynomial.
    !
    dx = xval - cs%x(j)
    yval = cs%y(j) + dx*(cs%b(j) + dx*(cs%c(j) + dx*cs%d(j)))
    
    return
  end subroutine lookupNonuniformSpline

  subroutine lookupUniformSpline(cs, xval, yval)
    
    implicit none

    type (cubicSpline), intent(in) :: cs
    real( kind = DP ), intent(in)  :: xval
    real( kind = DP ), intent(out) :: yval
    real( kind = DP ) ::  dx
    integer :: i, j
    !
    !  Find the interval J = [ cs%x(J), cs%x(J+1) ] that contains 
    !  or is nearest to xval.
    
    j = MAX(1, MIN(cs%n-1, int((xval-cs%x(1)) * cs%dx_i) + 1))
    
    dx = xval - cs%x(j)
    yval = cs%y(j) + dx*(cs%b(j) + dx*(cs%c(j) + dx*cs%d(j)))
    
    return
  end subroutine lookupUniformSpline

  subroutine lookupUniformSpline1d(cs, xval, yval, dydx)
    
    implicit none

    type (cubicSpline), intent(in) :: cs
    real( kind = DP ), intent(in)  :: xval
    real( kind = DP ), intent(out) :: yval, dydx
    real( kind = DP ) :: dx
    integer :: i, j
    
    !  Find the interval J = [ cs%x(J), cs%x(J+1) ] that contains 
    !  or is nearest to xval.


    j = MAX(1, MIN(cs%n-1, int((xval-cs%x(1)) * cs%dx_i) + 1))
    
    dx = xval - cs%x(j)
    yval = cs%y(j) + dx*(cs%b(j) + dx*(cs%c(j) + dx*cs%d(j)))

    dydx = cs%b(j) + dx*(2.0_dp * cs%c(j) + 3.0_dp * dx * cs%d(j))
       
    return
  end subroutine lookupUniformSpline1d

  subroutine lookupSpline(cs, xval, yval)

    type (cubicSpline), intent(in) :: cs
    real( kind = DP ), intent(inout) :: xval
    real( kind = DP ), intent(inout) :: yval
    
    if (cs%isUniform) then
       call lookupUniformSpline(cs, xval, yval)
    else
       call lookupNonuniformSpline(cs, xval, yval)
    endif

    return
  end subroutine lookupSpline
  
end module interpolation
