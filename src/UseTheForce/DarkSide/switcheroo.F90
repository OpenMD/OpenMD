!!
!! Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
!!$
module switcheroo

  use definitions
  use interpolation
  use status

  implicit none
  PRIVATE

#define __FORTRAN90
#include "UseTheForce/DarkSide/fSwitchingFunctionType.h"
  
  !! number of points for the spline approximations
  INTEGER, PARAMETER :: np = 150

  real ( kind = dp ), save :: rin
  real ( kind = dp ), save :: rout
  real ( kind = dp ), save :: rin2
  real ( kind = dp ), save :: rout2

  logical, save :: haveSplines = .false.
  logical, save :: switchIsCubic = .true.
  integer, save :: functionType = CUBIC

  ! spline structures
  type(cubicSpline), save :: r2Spline
  type(cubicSpline), save :: switchSpline

  public::set_switch_type
  public::set_switch
  public::get_switch
  public::delete_switch

contains

  subroutine set_switch(rinner, router)

    real ( kind = dp ), intent(in):: rinner, router
    real ( kind = dp ), dimension(np) :: xvals, yvals
    real ( kind = dp ), dimension(2) :: rCubVals, sCubVals
    real ( kind = dp ) :: rval, rval2, rval3, rval4, rval5
    real ( kind = dp ) :: rvaldi, rvaldi2, rvaldi3, rvaldi4, rvaldi5
    real ( kind = dp ) :: c0, c3, c4, c5, dx, r, r2
    integer :: i

    if (router .lt. rinner) then
       call handleError("set_switch", "router is less than rinner")
       return
    endif

    if ((router .lt. 0.0_dp) .or. (rinner .lt. 0.0_dp))  then
       call handleError("set_switch", "one of the switches is negative!")
       return
    endif

    rin = rinner
    rout = router
    rin2 = rinner * rinner
    rout2 = router * router

    if ((router-rinner) .lt. 1e-8)  then
       ! no reason to set up splines if the switching region is tiny
       return
    endif

    dx = (rout2-rin2) / dble(np-1)
    
    do i = 1, np
       r2 = rin2 + dble(i-1)*dx
       xvals(i) = r2
       yvals(i) = sqrt(r2)
    enddo

    call newSpline(r2spline, xvals, yvals, .true.)

    if (functionType .eq. FIFTH_ORDER_POLY) then
       c0 = 1.0_dp
       c3 = -10.0_dp
       c4 = 15.0_dp
       c5 = -6.0_dp

       dx = (rout-rin) / dble(np-1)
    
       do i = 1, np
          r = rin + dble(i-1)*dx
          xvals(i) = r

          rval = ( r - rin )
          rval2 = rval*rval
          rval3 = rval2*rval
          rval4 = rval2*rval2
          rval5 = rval3*rval2
          rvaldi = 1.0_dp/( rout - rin )
          rvaldi2 = rvaldi*rvaldi
          rvaldi3 = rvaldi2*rvaldi
          rvaldi4 = rvaldi2*rvaldi2
          rvaldi5 = rvaldi3*rvaldi2
          yvals(i)= c0 + c3*rval3*rvaldi3 + c4*rval4*rvaldi4 + c5*rval5*rvaldi5
       enddo
       
       call newSpline(switchSpline, xvals, yvals, .true.)
       
       switchIsCubic = .false.
    else
       rCubVals(1) = rin
       rCubVals(2) = rout
       sCubVals(1) = 1.0_dp
       sCubVals(2) = 0.0_dp      
       call newSpline(switchSpline, rCubVals, sCubVals, .true.)
    endif
    
    haveSplines = .true.
    return
  end subroutine set_switch

  subroutine set_switch_type(functionForm)
    integer, intent(in) :: functionForm    
    functionType = functionForm

    if ((functionType.eq.FIFTH_ORDER_POLY).or.(functionType.eq.CUBIC)) then
       if (haveSplines) then
          call delete_switch()
          call set_switch(rin, rout)
       endif
    else
       call handleError("set_switch_type", &
            "Unknown type of switching function!")
       return      
    endif
  end subroutine set_switch_type
  
  subroutine delete_switch()
    call deleteSpline(switchSpline)
    call deleteSpline(r2spline)
    return
  end subroutine delete_switch
  
  subroutine get_switch(r2, sw, dswdr, r, in_switching_region)

    real( kind = dp ), intent(in) :: r2
    real( kind = dp ), intent(inout) :: sw, dswdr, r
    logical, intent(inout) :: in_switching_region
    integer :: j
    real ( kind = dp ) :: a, b, c, d, dx

    sw = 1.0_dp
    dswdr = 0.0_dp
    in_switching_region = .false.

    if (r2.gt.rin2) then
       if (r2.gt.rout2) then

          sw = 0.0_dp
          dswdr = 0.0_dp
          return
          
       else         
          
          call lookupUniformSpline(r2Spline, r2, r)
          if (switchIsCubic) then
             ! super zippy automated use of precomputed spline coefficients
             dx = r - rin
             sw = switchSpline%y(1) + dx*(dx*(switchSpline%c(1) + &
                  dx*switchSpline%d(1)))
             dswdr = dx*(2.0_dp * switchSpline%c(1) + &
                  3.0_dp * dx * switchSpline%d(1))
          else
             call lookupUniformSpline1d(switchSpline, r, sw, dswdr)
          endif
          
          in_switching_region = .true.
          
          return          
       endif
    else
       return
    endif
    
  end subroutine get_switch

end module switcheroo
