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

module switcheroo

  use definitions
  use interpolation

  implicit none
  PRIVATE

#define __FORTRAN90
#include "UseTheForce/fSwitchingFunction.h"
#include "UseTheForce/DarkSide/fSwitchingFunctionType.h"

  real ( kind = dp ), dimension(NSWITCHTYPES) :: rin
  real ( kind = dp ), dimension(NSWITCHTYPES) :: rout
  real ( kind = dp ), dimension(NSWITCHTYPES) :: rin2
  real ( kind = dp ), dimension(NSWITCHTYPES) :: rout2
  real ( kind = dp ), save :: c0, c1, c2, c3, c4, c5

  logical, dimension(NSWITCHTYPES) :: isOK
  logical, save :: haveFunctionType = .false.
  logical, save :: haveSqrtSpline = .false.
  logical, save :: useSpline = .true.
  integer, save :: functionType = CUBIC


  ! spline variables
  type(cubicSpline), save :: scoef
  real ( kind = dp ), dimension(SPLINE_SEGMENTS) :: xValues
  real ( kind = dp ), dimension(SPLINE_SEGMENTS) :: yValues
  real ( kind = dp ), save :: dSqrt1, dSqrtN, range, dX
  real ( kind = dp ), save :: lowerBound
  logical, save :: uniformSpline = .true.

  public::set_switch
  public::set_function_type
  public::get_switch

contains

  subroutine set_switch(SwitchType, rinner, router)

    real ( kind = dp ), intent(in):: rinner, router
    integer, intent(in) :: SwitchType
    integer :: i

    if (SwitchType .gt. NSWITCHTYPES) then
       write(default_error, *) &
            'set_switch:  not that many switch types! '
       return
    endif

    isOK(SwitchType) = .false.

    if (router .lt. rinner) then
       write(default_error, *) &
            'set_switch:  router is less than rinner '
       return
    endif

    if ((router .lt. 0.0d0) .or. (rinner .lt. 0.0d0))  then
       write(default_error, *) &
            'set_switch:  one of the switches is negative!'
       return
    endif

    rin(SwitchType) = rinner
    rout(SwitchType) = router
    rin2(SwitchType) = rinner * rinner
    rout2(SwitchType) = router * router
    isOK(SwitchType) = .true.

    if (.not.haveSqrtSpline) then
       ! fill arrays for building the spline
       lowerBound = 1.0d0 ! the smallest value expected for r2
       range = rout2(SwitchType) - lowerBound
       dX = range / (SPLINE_SEGMENTS - 1)
       
       ! the spline is bracketed by lowerBound and rout2 endpoints
       xValues(1) = lowerBound
       yValues(1) = dsqrt(lowerBound)
       do i = 1, SPLINE_SEGMENTS-1
          xValues(i+1) = i * dX
          yValues(i+1) = dsqrt( i * dX )
       enddo
       
       ! set the endpoint derivatives
       dSqrt1 = 1 / ( 2.0d0 * dsqrt( xValues(1) ) )
       dSqrtN = 1 / ( 2.0d0 * dsqrt( xValues(SPLINE_SEGMENTS) ) )

       ! call newSpline to fill the coefficient array
       call newSpline(scoef, xValues, yValues, dSqrt1, dSqrtN, uniformSpline)
       
    endif
    
  end subroutine set_switch

  subroutine set_function_type(functionForm)
    integer, intent(in) :: functionForm    
    functionType = functionForm

    if (functionType .eq. FIFTH_ORDER_POLY) then
       c0 = 1.0d0
       c1 = 0.0d0
       c2 = 0.0d0
       c3 = -10.0d0
       c4 = 15.0d0
       c5 = -6.0d0
    endif
  end subroutine set_function_type

  subroutine get_switch(r2, sw, dswdr, r, SwitchType, in_switching_region)

    real( kind = dp ), intent(in) :: r2
    real( kind = dp ), intent(inout) :: sw, dswdr, r
    real( kind = dp ) :: ron, roff
    real( kind = dp ) :: rval, rval2, rval3, rval4, rval5
    real( kind = dp ) :: rvaldi, rvaldi2, rvaldi3, rvaldi4, rvaldi5
    integer, intent(in)    :: SwitchType
    logical, intent(inout) :: in_switching_region

    sw = 0.0d0
    dswdr = 0.0d0
    in_switching_region = .false.

    if (.not.isOK(SwitchType)) then
       write(default_error, *) &
            'get_switch:  this switching function has not been set up!'
       return
    endif

    if (r2.lt.rout2(SwitchType)) then
       if (r2.lt.rin2(SwitchType)) then

          sw = 1.0d0
          dswdr = 0.0d0
          return

       else
          if (useSpline) then
             call lookup_uniform_spline(scoef, r2, r)
          else
             r = dsqrt(r2)
          endif

          ron = rin(SwitchType)
          roff = rout(SwitchType)

          if (functionType .eq. FIFTH_ORDER_POLY) then
             rval = ( r - ron )
             rval2 = rval*rval
             rval3 = rval2*rval
             rval4 = rval2*rval2
             rval5 = rval3*rval2
             rvaldi = 1.0d0/( roff - ron )
             rvaldi2 = rvaldi*rvaldi
             rvaldi3 = rvaldi2*rvaldi
             rvaldi4 = rvaldi2*rvaldi2
             rvaldi5 = rvaldi3*rvaldi2
             sw = c0 + c1*rval*rvaldi + c2*rval2*rvaldi2 + c3*rval3*rvaldi3 &
                  + c4*rval4*rvaldi4 + c5*rval5*rvaldi5
             dswdr = c1*rvaldi + 2.0d0*c2*rval*rvaldi2 &
                  + 3.0d0*c3*rval2*rvaldi3 + 4.0d0*c4*rval3*rvaldi4 &
                  + 5.0d0*c5*rval4*rvaldi5

          else
             sw = (roff + 2.0d0*r - 3.0d0*ron)*(roff-r)**2/ ((roff-ron)**3)
             dswdr = 6.0d0*(r*r - r*ron - r*roff +roff*ron)/((roff-ron)**3)

          endif
          in_switching_region = .true.
          return          
       endif
    else
       return
    endif

  end subroutine get_switch

end module switcheroo
