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

  implicit none
  PRIVATE

#define __FORTRAN90
#include "UseTheForce/fSwitchingFunction.h"

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
