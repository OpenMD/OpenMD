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
  is_Electrostatic = (atp%is_Charge .ne. 0) .or. (atp%is_Dipole .ne. 0) .and. &
                     (atp%is_Quadrupole .ne. 0)
  
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
