#ifdef __OPENMD_C
#ifndef __FMNMINTERACTIONS
#define __FMNMINTERACTIONS

#define MNM_NUM_MNM_TYPES 4
#define MNM_LENNARDJONES  1
#define MNM_REPULSIVEMORSE 2
#define MNM_SHIFTEDMORSE 3
#define MNM_MAW 4

// dummy integer below is a workaround for
// a bug in the intel compiler, which can't
// seem to pack structs the same way on C and
// fortran.

typedef  struct{
  int MNMInteractionType;
  int dummy;
  int metal_atid;
  int nonmetal_atid;
  RealType R0;
  RealType D0;
  RealType beta0;  
  RealType betaH;  
  RealType ca1;  
  RealType cb1;  
  RealType sigma;
  RealType epsilon;
} MNMtype;


#endif
#endif /*__OPENMD_C*/

#ifdef  __FORTRAN90

  INTEGER, PARAMETER:: MNM_NUM_MNM_TYPES = 4
  INTEGER, PARAMETER:: MNM_LENNARDJONES  = 1
  INTEGER, PARAMETER:: MNM_REPULSIVEMORSE = 2
  INTEGER, PARAMETER:: MNM_SHIFTEDMORSE = 3
  INTEGER, PARAMETER:: MNM_MAW = 4

  ! dummy integer below is a workaround for
  ! a bug in the intel compiler, which does not 
  ! seem to pack structs the same way on C and
  ! fortran.

  type :: MNMtype
    SEQUENCE
    integer :: MNMInteractionType
    integer :: dummy
    integer :: metal_atid
    integer :: nonmetal_atid
    real(kind=dp) :: R0
    real(kind=dp) :: D0
    real(kind=dp) :: beta0
    real(kind=dp) :: betaH
    real(kind=dp) :: ca1
    real(kind=dp) :: cb1
    real(kind=dp) :: sigma
    real(kind=dp) :: epsilon
  end type MNMtype

#endif
