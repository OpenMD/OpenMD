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

!! module defines atypes available to simulation

module atype_module

  use vector_class
  implicit none
  private

  type (Vector), pointer, public :: atypes => null()


  public :: new_atype
  public :: delete_atypes

  

contains

  subroutine new_atype(ident, is_Directional, is_LennardJones, &
       is_Electrostatic, is_Charge, is_Dipole, is_Quadrupole, &
       is_Sticky, is_StickyPower, is_GayBerne, is_EAM, is_Shape, &
       is_FLARB, is_SC, status)
    integer :: myATID, c_ident
    integer,intent(in) :: ident
    logical,intent(in) :: is_Directional, is_LennardJones, is_Electrostatic
    logical,intent(in) :: is_Charge, is_Dipole, is_Quadrupole
    logical,intent(in) :: is_Sticky, is_StickyPower, is_GayBerne, is_EAM
    logical,intent(in) :: is_Shape, is_FLARB, is_SC
    integer,intent(out) :: status

    integer :: me

    status = 0

    if (.not. associated(atypes)) then
       !! There are 17 properties to worry about for now.  
       !! Fix this if needed for more atomic properties
       atypes => initialize(18)
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
    call setElementProperty(atypes, me, "is_Quadrupole", is_Quadrupole)
    call setElementProperty(atypes, me, "is_Sticky", is_Sticky)
    call setElementProperty(atypes, me, "is_StickyPower", is_StickyPower)
    call setElementProperty(atypes, me, "is_GayBerne", is_GayBerne)
    call setElementProperty(atypes, me, "is_EAM", is_EAM)
    call setElementProperty(atypes, me, "is_Shape", is_Shape)
    call setElementProperty(atypes, me, "is_FLARB", is_FLARB)
    call setElementProperty(atypes, me, "is_SC", is_SC)

  end subroutine new_atype

  subroutine delete_atypes()
    atypes => destroy(atypes)
  end subroutine delete_atypes




end module atype_module

