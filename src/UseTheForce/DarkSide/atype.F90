!! module defines atypes available to simulation

module atype_module

  use vector_class
  implicit none
  private

  type (Vector), pointer, public :: atypes => null()

  public :: new_atype
  
contains
  
  subroutine new_atype(ident, is_Directional, is_LennardJones, &
       is_Electrostatic, is_Charge, is_Dipole, is_Sticky, is_GayBerne, &
       is_EAM, is_Shape, is_FLARB, status)
    
    integer,intent(in) :: ident
    logical,intent(in) :: is_Directional, is_LennardJones, is_Electrostatic
    logical,intent(in) :: is_Charge, is_Dipole, is_Sticky, is_GayBerne, is_EAM
    logical,intent(in) :: is_Shape, is_FLARB
    integer,intent(out) :: status

    integer :: me
        
    status = 0
    
    if (.not. associated(atypes)) then
       !! There are 16 properties to worry about for now.  
       !! Fix this if needed for more atomic properties
       atypes => initialize(17)
       if (.not.associated(atypes)) then
          status = -1
          return
       endif
    endif
    
    me = addElement(atypes)
    call setElementProperty(atypes, me, "c_ident", ident)
    
    call setElementProperty(atypes, me, "is_Directional", is_Directional)
    call setElementProperty(atypes, me, "is_LennardJones", is_LennardJones)
    call setElementProperty(atypes, me, "is_Electrostatic", is_Electrostatic)
    call setElementProperty(atypes, me, "is_Charge", is_Charge)
    call setElementProperty(atypes, me, "is_Dipole", is_Dipole)
    call setElementProperty(atypes, me, "is_Sticky", is_Sticky)
    call setElementProperty(atypes, me, "is_GayBerne", is_GayBerne)
    call setElementProperty(atypes, me, "is_EAM", is_EAM)
    call setElementProperty(atypes, me, "is_Shape", is_Shape)
    call setElementProperty(atypes, me, "is_FLARB", is_FLARB)
    
  end subroutine new_atype
  
end module atype_module

! provide interface for c calls....

subroutine makeatype(atp, status)
  
  use atype_module, ONLY: new_atype

#define __FORTRAN90
#include "types/AtomTypeProperties.h"    

  type(AtomTypeProperties), intent(in) :: atp
  integer, intent(inout) :: status

  integer :: ident
  logical :: is_Directional, is_LennardJones, is_Electrostatic
  logical :: is_Charge, is_Dipole, is_Sticky, is_GayBerne, is_EAM
  logical :: is_Shape, is_FLARB
  
  ident = atp%ident
  is_Directional = (atp%is_Directional .ne. 0)
  is_LennardJones = (atp%is_LennardJones .ne. 0)
  is_Electrostatic = (atp%is_Electrostatic .ne. 0)
  is_Charge = (atp%is_Charge .ne. 0)
  is_Dipole = (atp%is_Dipole .ne. 0)
  is_Sticky = (atp%is_Sticky .ne. 0)
  is_GayBerne = (atp%is_GayBerne .ne. 0)
  is_EAM = (atp%is_EAM .ne. 0)
  is_Shape = (atp%is_Shape .ne. 0)
  is_FLARB = (atp%is_FLARB .ne. 0)

  call new_atype(ident, is_Directional, is_LennardJones, is_Electrostatic, &
       is_Charge, is_Dipole, is_Sticky, is_GayBerne, is_EAM, is_Shape, &
       is_FLARB, status)
  
end subroutine
