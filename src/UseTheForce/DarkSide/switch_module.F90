module switcheroo

  use definitions

  implicit none
  PRIVATE

#define __FORTRAN90
#include "fSwitchingFunction.h"

  real ( kind = dp ), dimension(NSWITCHTYPES) :: rin
  real ( kind = dp ), dimension(NSWITCHTYPES) :: rout
  real ( kind = dp ), dimension(NSWITCHTYPES) :: rin2
  real ( kind = dp ), dimension(NSWITCHTYPES) :: rout2

  logical, dimension(NSWITCHTYPES) :: isOK


  public::set_switch
  public::get_switch

contains

  subroutine set_switch(SwitchType, rinner, router)

    real ( kind = dp ), intent(in):: rinner, router
    integer, intent(in) :: SwitchType

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

  end subroutine set_switch

  subroutine get_switch(r2, sw, dswdr, r, SwitchType, in_switching_region)

    real( kind = dp ), intent(in) :: r2
    real( kind = dp ), intent(inout) :: sw, dswdr, r
    real( kind = dp ) :: ron, roff
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
          
          r = dsqrt(r2)
          
          ron = rin(SwitchType)
          roff = rout(SwitchType)
          
          sw = (roff + 2.0d0*r - 3.0d0*ron)*(roff-r)**2/ ((roff-ron)**3)
          dswdr = 6.0d0*(r*r - r*ron - r*roff +roff*ron)/((roff-ron)**3)
          
          in_switching_region = .true.
          return          
       endif
    else
       return
    endif    
       
  end subroutine get_switch
end module switcheroo
