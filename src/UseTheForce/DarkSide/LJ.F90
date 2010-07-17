module ljmod
  use, intrinsic :: iso_c_binding
  use definitions
  implicit none
  private
  type(c_funptr) :: ljClass
  save

  public :: getSigma
  public :: getEpsilon
  public :: setLJDefaultCutoff
  public :: do_lj_pair

  real(kind=dp) function getSigma(atid)
    integer :: atid
    re


  function LJgetSigma(ljClass, atid) result(sigma) bind(C,NAME="LJ_getSigma")
    use iso_c_binding
    type(c_funptr), value :: ljClass
    integer(c_int), value :: atid
    real(c_double) :: sigma
  end function getSigma

  function LJgetEpsilon(ljClass, atid) result(epsilon) bind(C,NAME="LJ_getEpsilon")
    use iso_c_binding
    type(c_funptr), value :: ljClass
    integer(c_int), value :: atid
    real(c_double) :: epsilon
  end function getEpsilon

    
end module ljmod
     
