!! module defines atypes available to simulation

module atype_module
  use definitions, only: dp
  use vector_class
  implicit none
  private
  
  type (Vector), pointer, public :: atypes => null()
  
  
  public :: new_atype
  
contains
  
  subroutine new_atype(c_ident, is_LJ, is_Sticky, is_DP, is_GB, &
       is_EAM, is_Charge, lj_epsilon, lj_sigma, charge, dipole_moment, &
       status)
    
    real( kind = dp ), intent(in) :: lj_epsilon
    real( kind = dp ), intent(in) :: lj_sigma
    real( kind = dp ), intent(in) :: dipole_moment
    real( kind = dp ), intent(in) :: charge

    integer, intent(in)  :: c_ident
    integer, intent(out) :: status
    integer, intent(in)  :: is_Sticky
    integer, intent(in)  :: is_DP
    integer, intent(in)  :: is_GB
    integer, intent(in)  :: is_EAM
    integer, intent(in)  :: is_LJ
    integer, intent(in)  :: is_Charge
    integer :: me
    logical :: l_is_LJ, l_is_DP, l_is_Sticky, l_is_GB
    logical :: l_is_EAM, l_is_Charge
    integer :: FFcheckStatus
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
    call setElementProperty(atypes, me, "c_ident", c_ident)

    l_is_LJ = (is_LJ .ne. 0)
    l_is_DP = (is_DP .ne. 0)
    l_is_Sticky = (is_Sticky .ne. 0)
    l_is_GB = (is_GB .ne. 0)
    l_is_EAM = (is_EAM .ne. 0)
    l_is_Charge = (is_Charge .ne. 0)

    call setElementProperty(atypes, me, "is_LJ", l_is_LJ)
    call setElementProperty(atypes, me, "is_DP", l_is_DP)
    call setElementProperty(atypes, me, "is_Sticky", l_is_Sticky)
    call setElementProperty(atypes, me, "is_GB", l_is_GB)
    call setElementProperty(atypes, me, "is_EAM", l_is_EAM)
    call setElementProperty(atypes, me, "is_Charge", l_is_Charge)

    if (l_is_LJ) then
       call setElementProperty(atypes, me, "lj_sigma", lj_sigma)
       call setElementProperty(atypes, me, "lj_epsilon", lj_epsilon)
    endif
    if (l_is_DP) then
       call setElementProperty(atypes, me, "dipole_moment", dipole_moment)
    endif
    if (l_is_Charge) then
       call setElementProperty(atypes, me, "charge", charge)
    endif

  end subroutine new_atype

end module atype_module
  ! provide interface for c calls....
subroutine makeatype(c_ident, is_LJ, is_Sticky, is_DP, is_GB, &
       is_EAM, is_Charge, lj_epsilon, lj_sigma, charge, dipole_moment, &
       status)
    use definitions, only: dp   
    use atype_module, ONLY: new_atype
    
    real( kind = dp ), intent(in) :: lj_epsilon
    real( kind = dp ), intent(in) :: lj_sigma
    real( kind = dp ), intent(in) :: dipole_moment
    real( kind = dp ), intent(in) :: charge

    integer, intent(in)  :: c_ident
    integer, intent(out) :: status
    integer, intent(in)  :: is_Sticky
    integer, intent(in)  :: is_DP
    integer, intent(in)  :: is_GB
    integer, intent(in)  :: is_EAM
    integer, intent(in)  :: is_LJ
    integer, intent(in)  :: is_Charge

    call module_new_atype(c_ident, is_LJ, is_Sticky, is_DP, is_GB, &
       is_EAM, is_Charge, lj_epsilon, lj_sigma, charge, dipole_moment, &
       status)
       
end subroutine