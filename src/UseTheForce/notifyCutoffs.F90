module notifyCutoffs
  
  use definitions
  use do_Forces, only:      setRlistDF
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
      
      call setRlistDF( rlist )
      call setCutoffsRF( rcut, rsw )
      call setCutoffLJ( rcut, do_shift, localError )
      call setCutoffEAM( rcut, localError)
      call set_switch(GROUP_SWITCH, rsw, rcut)

    end subroutine cutoffNotify

  end module notifyCutoffs
  
  subroutine notifyFortranCutoffs(this_rcut, this_rsw, this_rlist )
    use notifyCutoffs
    use definitions, ONLY : dp
    
    real(kind=dp), intent(in) :: this_rcut, this_rsw, this_rlist
    
    call cutoffNotify(this_rcut, this_rsw, this_rlist )
    
  end subroutine 
