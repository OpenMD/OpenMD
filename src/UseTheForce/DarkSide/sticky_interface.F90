subroutine newStickyType(c_ident, w0, v0, v0p, rl, ru, rlp, rup, isError)
  
  use definitions, ONLY : dp   
  use sticky, ONLY : module_newStickyType => newStickyType
  
  integer, intent(inout) :: c_ident, isError
  real( kind = dp ), intent(inout) :: w0, v0, v0p, rl, ru, rlp, rup
  
  call module_newStickyType(c_ident, w0, v0, v0p, rl, ru, rlp, rup, &
       isError)
  
end subroutine newStickyType

subroutine destroyStickyTypes()
  use sticky, ONLY : module_destroyStickyTypes => destroyStickyTypes
  call module_destroyStickyTypes()

end subroutine destroyStickyTypes
