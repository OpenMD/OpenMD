! provide interface for c calls....

subroutine makeatype(atp, status)
  
  use atype_module, ONLY: new_atype

#define __FORTRAN90
#include "types/AtomTypeProperties.h"    

  type(AtomTypeProperties), intent(in) :: atp
  integer, intent(inout) :: status

  integer :: ident
  logical :: is_Directional, is_LennardJones, is_Electrostatic
  logical :: is_Charge, is_Dipole, is_Quadrupole
  logical :: is_Sticky, is_GayBerne, is_EAM, is_Shape, is_FLARB

  ident = atp%ident
  is_Directional = (atp%is_Directional .ne. 0)
  is_LennardJones = (atp%is_LennardJones .ne. 0)
  is_Electrostatic = ((atp%is_Charge .ne. 0) .or. (atp%is_Dipole .ne. 0)) &
       .or. (atp%is_Quadrupole .ne. 0)
  
  is_Charge = (atp%is_Charge .ne. 0)
  is_Dipole = (atp%is_Dipole .ne. 0)
  is_Quadrupole = (atp%is_Quadrupole .ne. 0)
  is_Sticky = (atp%is_Sticky .ne. 0)
  is_GayBerne = (atp%is_GayBerne .ne. 0)
  is_EAM = (atp%is_EAM .ne. 0)
  is_Shape = (atp%is_Shape .ne. 0)
  is_FLARB = (atp%is_FLARB .ne. 0)

  call new_atype(ident, is_Directional, is_LennardJones, is_Electrostatic, &
       is_Charge, is_Dipole, is_Quadrupole, is_Sticky, is_GayBerne, is_EAM, &
       is_Shape, is_FLARB, status)
  
end subroutine

! clears memory up
subroutine deleteAtypes()
  use atype_module,ONLY: delete_atypes
  call delete_atypes()
end subroutine deleteatypes
