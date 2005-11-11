subroutine setFunctionType(the_FT)
  use switcheroo, ONLY : module_setFT => set_function_type
  integer,intent(inout) :: the_FT
  call module_setFT(the_FT)
end subroutine setFunctionType

