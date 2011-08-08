! provide interface for c calls....

subroutine addMNMInteraction(myInteraction)

  use MetalNonMetal, ONLY: module_addInteraction => addInteraction
  use definitions, ONLY : dp
#define __FORTRAN90
#include "UseTheForce/DarkSide/fMnMInteractions.h"    

  type(MNMType), intent(inout) :: myInteraction
  call module_addInteraction(myInteraction)

end subroutine addMNMInteraction

! clears memory up
subroutine deleteMNMInteractions()
  use MetalNonMetal,ONLY: module_deleteInteractions => deleteInteractions
  call module_deleteInteractions()
end subroutine deleteMNMInteractions
