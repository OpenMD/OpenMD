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
!! @version $Id: interpolation.F90,v 1.1 2006-04-14 19:57:04 gezelter Exp $


module  INTERPOLATION
  use definitions
  use status
  implicit none
  PRIVATE

  character(len = statusMsgSize) :: errMSG

  type, public :: cubicSpline
     private
     integer :: np = 0
     real(kind=dp) :: dx
     real(kind=dp) :: dx_i
     real (kind=dp), pointer,dimension(:)   :: x => null()
     real (kind=dp), pointer,dimension(4,:) :: c => null()
  end type cubicSpline

  interface splineLookup
     module procedure multiSplint
     module procedure splintd
     module procedure splintd1
     module procedure splintd2
  end interface

  interface newSpline
     module procedure newSplineWithoutDerivs
     module procedure newSplineWithDerivs
  end interface

  public :: deleteSpline

contains


  subroutine newSplineWithoutDerivs(cs, x, y, yp1, ypn, boundary)

    !************************************************************************
    !
    ! newSplineWithoutDerivs solves for slopes defining a cubic spline.
    !
    !  Discussion:
    !
    !    A tridiagonal linear system for the unknown slopes S(I) of
    !    F at x(I), I=1,..., N, is generated and then solved by Gauss
    !    elimination, with S(I) ending up in cs%C(2,I), for all I.
    !
    !  Reference:
    !
    !    Carl DeBoor,
    !    A Practical Guide to Splines,
    !    Springer Verlag.
    !
    !  Parameters:
    !
    !    Input, real x(N), the abscissas or X values of
    !    the data points.  The entries of TAU are assumed to be
    !    strictly increasing.
    !
    !    Input, real y(I), contains the function value at x(I) for 
    !      I = 1, N.
    !
    !    yp1 contains the slope at x(1) and ypn contains
    !    the slope at x(N).
    !
    !    On output, the intermediate slopes at x(I) have been
    !    stored in cs%C(2,I), for I = 2 to N-1.

    implicit none

    type (cubicSpline), intent(inout) :: cs
    real( kind = DP ), intent(in) :: x(:), y(:)
    real( kind = DP ), intent(in) :: yp1, ypn
    character(len=*), intent(in) :: boundary
    real( kind = DP ) :: g, divdif1, divdif3, dx
    integer :: i, alloc_error, np

    alloc_error = 0

    if (cs%np .ne. 0) then
       call handleWarning("interpolation::newSplineWithoutDerivs", &
            "Type was already created")
       call deleteSpline(cs)
    end if

    ! make sure the sizes match

    if (size(x) .ne. size(y)) then
       call handleError("interpolation::newSplineWithoutDerivs", &
            "Array size mismatch")
    end if

    np = size(x)
    cs%np = np

    allocate(cs%x(np), stat=alloc_error)
    if(alloc_error .ne. 0) then
       call handleError("interpolation::newSplineWithoutDerivs", &
            "Error in allocating storage for x")
    endif

    allocate(cs%c(4,np), stat=alloc_error)
    if(alloc_error .ne. 0) then
       call handleError("interpolation::newSplineWithoutDerivs", &
            "Error in allocating storage for c")
    endif
       
    do i = 1, np
       cs%x(i) = x(i)
       cs%c(1,i) = y(i)       
    enddo

    if ((boundary.eq.'l').or.(boundary.eq.'L').or. &
         (boundary.eq.'b').or.(boundary.eq.'B')) then
       cs%c(2,1) = yp1
    else
       cs%c(2,1) = 0.0_DP
    endif
    if ((boundary.eq.'u').or.(boundary.eq.'U').or. &
         (boundary.eq.'b').or.(boundary.eq.'B')) then
       cs%c(2,1) = ypn
    else
       cs%c(2,1) = 0.0_DP
    endif

    !
    !  Set up the right hand side of the linear system.
    !
    do i = 2, cs%np - 1
       cs%c(2,i) = 3.0_DP * ( &
            (x(i) - x(i-1)) * (cs%c(1,i+1) - cs%c(1,i)) / (x(i+1) - x(i)) + &
            (x(i+1) - x(i)) * (cs%c(1,i) - cs%c(1,i-1)) / (x(i) - x(i-1)))
    end do
    !
    !  Set the diagonal coefficients.
    !
    cs%c(4,1) = 1.0_DP
    do i = 2, cs%np - 1 
       cs%c(4,i) = 2.0_DP * ( x(i+1) - x(i-1) )
    end do
    cs%c(4,n) = 1.0_DP
    !
    !  Set the off-diagonal coefficients.
    !
    cs%c(3,1) = 0.0_DP
    do i = 2, cs%np
       cs%c(3,i) = x(i) - x(i-1)
    end do
    !
    !  Forward elimination.
    !
    do i = 2, cs%np - 1
       g = -cs%c(3,i+1) / cs%c(4,i-1)
       cs%c(4,i) = cs%c(4,i) + g * cs%c(3,i-1)
       cs%c(2,i) = cs%c(2,i) + g * cs%c(2,i-1)
    end do
    !
    !  Back substitution for the interior slopes.
    !
    do i = cs%np - 1, 2, -1
       cs%c(2,i) = ( cs%c(2,i) - cs%c(3,i) * cs%c(2,i+1) ) / cs%c(4,i)
    end do
    !
    !  Now compute the quadratic and cubic coefficients used in the 
    !  piecewise polynomial representation.
    !
    do i = 1, cs%np - 1
       dx = x(i+1) - x(i)
       divdif1 = ( cs%c(1,i+1) - cs%c(1,i) ) / dx
       divdif3 = cs%c(2,i) + cs%c(2,i+1) - 2.0_DP * divdif1
       cs%c(3,i) = ( divdif1 - cs%c(2,i) - divdif3 ) / dx
       cs%c(4,i) = divdif3 / ( dx * dx )
    end do

    cs%c(3,np) = 0.0_DP
    cs%c(4,np) = 0.0_DP

    cs%dx = dx
    cs%dxi = 1.0_DP / dx
    return
  end subroutine newSplineWithoutDerivs

  subroutine newSplineWithDerivs(cs, x, y, yp)

    !************************************************************************
    !
    ! newSplineWithDerivs 

    implicit none

    type (cubicSpline), intent(inout) :: cs
    real( kind = DP ), intent(in) :: x(:), y(:), yp(:)
    real( kind = DP ) :: g, divdif1, divdif3, dx
    integer :: i, alloc_error, np

    alloc_error = 0

    if (cs%np .ne. 0) then
       call handleWarning("interpolation::newSplineWithDerivs", &
            "Type was already created")
       call deleteSpline(cs)
    end if

    ! make sure the sizes match

    if ((size(x) .ne. size(y)).or.(size(x) .ne. size(yp))) then
       call handleError("interpolation::newSplineWithDerivs", &
            "Array size mismatch")
    end if
    
    np = size(x)
    cs%np = np

    allocate(cs%x(np), stat=alloc_error)
    if(alloc_error .ne. 0) then
       call handleError("interpolation::newSplineWithDerivs", &
            "Error in allocating storage for x")
    endif
    
    allocate(cs%c(4,np), stat=alloc_error)
    if(alloc_error .ne. 0) then
       call handleError("interpolation::newSplineWithDerivs", &
            "Error in allocating storage for c")
    endif
    
    do i = 1, np
       cs%x(i) = x(i)
       cs%c(1,i) = y(i)       
       cs%c(2,i) = yp(i)
    enddo
    !
    !  Set the diagonal coefficients.
    !
    cs%c(4,1) = 1.0_DP
    do i = 2, cs%np - 1 
       cs%c(4,i) = 2.0_DP * ( x(i+1) - x(i-1) )
    end do
    cs%c(4,n) = 1.0_DP
    !
    !  Set the off-diagonal coefficients.
    !
    cs%c(3,1) = 0.0_DP
    do i = 2, cs%np
       cs%c(3,i) = x(i) - x(i-1)
    end do
    !
    !  Forward elimination.
    !
    do i = 2, cs%np - 1
       g = -cs%c(3,i+1) / cs%c(4,i-1)
       cs%c(4,i) = cs%c(4,i) + g * cs%c(3,i-1)
       cs%c(2,i) = cs%c(2,i) + g * cs%c(2,i-1)
    end do
    !
    !  Back substitution for the interior slopes.
    !
    do i = cs%np - 1, 2, -1
       cs%c(2,i) = ( cs%c(2,i) - cs%c(3,i) * cs%c(2,i+1) ) / cs%c(4,i)
    end do
    !
    !  Now compute the quadratic and cubic coefficients used in the 
    !  piecewise polynomial representation.
    !
    do i = 1, cs%np - 1
       dx = x(i+1) - x(i)
       divdif1 = ( cs%c(1,i+1) - cs%c(1,i) ) / dx
       divdif3 = cs%c(2,i) + cs%c(2,i+1) - 2.0_DP * divdif1
       cs%c(3,i) = ( divdif1 - cs%c(2,i) - divdif3 ) / dx
       cs%c(4,i) = divdif3 / ( dx * dx )
    end do

    cs%c(3,np) = 0.0_DP
    cs%c(4,np) = 0.0_DP

    cs%dx = dx
    cs%dxi = 1.0_DP / dx

    return
  end subroutine newSplineWithoutDerivs

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
    
    this%np = 0
    
  end subroutine deleteSpline

  subroutine lookup_nonuniform_spline(cs, xval, yval)
    
    !*************************************************************************
    !
    ! lookup_nonuniform_spline evaluates a piecewise cubic Hermite interpolant.
    !
    !  Discussion:
    !
    !    newSpline must be called first, to set up the
    !    spline data from the raw function and derivative data.
    !
    !  Modified:
    !
    !    06 April 1999
    !
    !  Reference:
    !
    !    Conte and de Boor,
    !    Algorithm PCUBIC,
    !    Elementary Numerical Analysis, 
    !    1973, page 234.
    !
    !  Parameters:
    !
    implicit none

    type (cubicSpline), intent(in) :: cs
    real( kind = DP ), intent(in)  :: xval
    real( kind = DP ), intent(out) :: yval
    integer :: i, j
    !
    !  Find the interval J = [ cs%x(J), cs%x(J+1) ] that contains 
    !  or is nearest to xval.
    !
    j = cs%np - 1

    do i = 1, cs%np - 2

       if ( xval < cs%x(i+1) ) then
          j = i
          exit
       end if

    end do
    !
    !  Evaluate the cubic polynomial.
    !
    dx = xval - cs%x(j)

    yval = cs%c(1,j) + dx * ( cs%c(2,j) + dx * ( cs%c(3,j) + dx * cs%c(4,j) ) )
    
    return
  end subroutine lookup_nonuniform_spline

  subroutine lookup_uniform_spline(cs, xval, yval)
    
    !*************************************************************************
    !
    ! lookup_uniform_spline evaluates a piecewise cubic Hermite interpolant.
    !
    !  Discussion:
    !
    !    newSpline must be called first, to set up the
    !    spline data from the raw function and derivative data.
    !
    !  Modified:
    !
    !    06 April 1999
    !
    !  Reference:
    !
    !    Conte and de Boor,
    !    Algorithm PCUBIC,
    !    Elementary Numerical Analysis, 
    !    1973, page 234.
    !
    !  Parameters:
    !
    implicit none

    type (cubicSpline), intent(in) :: cs
    real( kind = DP ), intent(in)  :: xval
    real( kind = DP ), intent(out) :: yval
    integer :: i, j
    !
    !  Find the interval J = [ cs%x(J), cs%x(J+1) ] that contains 
    !  or is nearest to xval.

    j = MAX(1, MIN(cs%np, idint((xval-cs%x(1)) * cs%dxi) + 1))

    dx = xval - cs%x(j)

    yval = cs%c(1,j) + dx * ( cs%c(2,j) + dx * ( cs%c(3,j) + dx * cs%c(4,j) ) )
    
    return
  end subroutine lookup_uniform_spline
  
end module INTERPOLATION
