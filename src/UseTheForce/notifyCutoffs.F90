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

module notifyCutoffs

  use definitions
  use doForces, only:       createRCuts
  use reaction_field, only: setCutoffsRF
  use lj, only:             setCutoffLJ
  use eam, only:            setCutoffEAM
  use switcheroo, only:     set_switch
  use status
  implicit none

  PRIVATE

  character(len = statusMsgSize) :: errMsg

#define __FORTRAN90
#include "UseTheForce/fSwitchingFunction.h"

  public::cutoffNotify

contains

  subroutine cutoffNotify( this_rcut, this_rsw, this_rlist )

    real(kind=dp), intent(in) :: this_rcut, this_rsw, this_rlist

    real(kind=dp) :: rsw, rcut, rlist
    integer :: localError
    logical :: do_shift

    rcut   = this_rcut
    rsw    = this_rsw
    rlist  = this_rlist

    if (rcut .lt. rsw) then

       write(errMsg, *) 'cutoffRadius is ', rcut, newline // tab, &            
            'but switchingRadius is set larger at ', rsw , newline // tab, &
            'That is probably not what you wanted to do!'

       call handleWarning("cutoffNotify", errMsg)

    endif

    if (rlist .lt. rcut) then

       write(errMsg, *) 'neighborListRadius is ', rlist, newline &
            // tab,  'but cutoffRadius is set larger at ', rcut , newline &
            // tab,  'That is probably a programming error!'

       call handleWarning("cutoffNotify", errMsg)

    endif

    do_shift = .false.
    if (abs(rcut-rsw) .lt. 0.0001) then

       write(errMsg, *) &
            'cutoffRadius and switchingRadius are set to the same', newline &
            // tab, 'value.  OOPSE will use shifted Lennard-Jones', newline &
            // tab, 'potentials instead of switching functions.'

       call handleInfo("cutoffNotify", errMsg)

       do_shift = .true.

    endif

    call createRCuts( defaultRlist=rlist,stat=localError )
    call setCutoffsRF( rcut, rsw )
    call setCutoffLJ( rcut, do_shift, localError )
    call setCutoffEAM( rcut, localError)
    call set_switch(GROUP_SWITCH, rsw, rcut)

    if (localError /= 0) then
       write(errMsg, *) 'An error has occured in setting the default cutoff'
       call handleError("cutoffNotify", errMsg)
    end if


  end subroutine cutoffNotify

end module notifyCutoffs
