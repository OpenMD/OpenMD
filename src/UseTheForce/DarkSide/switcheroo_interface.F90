subroutine setFunctionType(the_FT)
  use switcheroo, ONLY : module_setFT => set_switch_type
  integer,intent(inout) :: the_FT
  call module_setFT(the_FT)
end subroutine setFunctionType

subroutine deleteSwitch()
  use switcheroo, ONLY : module_delete_switch => delete_switch
  call module_delete_switch()
end subroutine deleteSwitch

