subroutine newChargeType(ident, charge, status)

  use charge_charge, ONLY : module_newChargeType => newChargeType

  integer, parameter :: DP = selected_real_kind(15)
  integer,intent(inout) :: ident
  real(kind=dp),intent(inout) :: charge
  integer,intent(inout) :: status
  
  call module_newChargeType(ident, charge, status)
  
end subroutine newChargeType
