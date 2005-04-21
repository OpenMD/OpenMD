subroutine makeShape(nContactFuncs, ContactFuncLValue, &
     ContactFuncMValue, ContactFunctionType, ContactFuncCoefficient, &
     nRangeFuncs, RangeFuncLValue, RangeFuncMValue, RangeFunctionType, &
     RangeFuncCoefficient, nStrengthFuncs, StrengthFuncLValue, &
     StrengthFuncMValue, StrengthFunctionType, StrengthFuncCoefficient, &
     myATID, status)

  use definitions
  use shapes, only: newShapeType

  integer :: nContactFuncs 
  integer :: nRangeFuncs 
  integer :: nStrengthFuncs 
  integer :: status
  integer :: c_ident

  integer, dimension(nContactFuncs) :: ContactFuncLValue           
  integer, dimension(nContactFuncs) :: ContactFuncMValue           
  integer, dimension(nContactFuncs) :: ContactFunctionType         
  real(kind=dp), dimension(nContactFuncs) :: ContactFuncCoefficient
  integer, dimension(nRangeFuncs) :: RangeFuncLValue             
  integer, dimension(nRangeFuncs) :: RangeFuncMValue             
  integer, dimension(nRangeFuncs) :: RangeFunctionType           
  real(kind=dp), dimension(nRangeFuncs) :: RangeFuncCoefficient  
  integer, dimension(nStrengthFuncs) :: StrengthFuncLValue          
  integer, dimension(nStrengthFuncs) :: StrengthFuncMValue          
  integer, dimension(nStrengthFuncs) :: StrengthFunctionType        
  real(kind=dp), dimension(nStrengthFuncs) :: StrengthFuncCoefficient

  call newShapeType(nContactFuncs, ContactFuncLValue, &
       ContactFuncMValue, ContactFunctionType, ContactFuncCoefficient, &
       nRangeFuncs, RangeFuncLValue, RangeFuncMValue, RangeFunctionType, &
       RangeFuncCoefficient, nStrengthFuncs, StrengthFuncLValue, &
       StrengthFuncMValue, StrengthFunctionType, StrengthFuncCoefficient, &
       c_ident, status)

  return
end subroutine makeShape

subroutine completeShapeFF(status)

  use shapes, only: complete_Shape_FF

  integer, intent(out)  :: status
  integer :: myStatus

  myStatus = 0

  call complete_Shape_FF(myStatus)

  status = myStatus

  return
end subroutine completeShapeFF

subroutine destroyShapeTypes()
  use shapes, only: module_destroyShapeTypes => destroyShapeTypes
  call module_destroyShapeTypes()

end subroutine destroyShapeTypes
