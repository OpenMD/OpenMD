#ifdef __C
#ifndef __FSIMULATION
#define __FSIMULATION
/** This header provides dual access for the simulation structure between 
    fortran and C for the simtype structure. NOTE: Sequence of struct 
    components must match between C and fortran and in general be packed 
    double,int,char. 
*/
typedef  struct{
  double dielect;
  int SIM_uses_PBC;
  int SIM_uses_LJ;
  int SIM_uses_sticky;
  int SIM_uses_charges;
  int SIM_uses_dipoles;
  int SIM_uses_RF;
  int SIM_uses_GB;
  int SIM_uses_EAM;
} simtype;
#endif //__FSIMULATION
#endif //__C

#ifdef  __FORTRAN90

type, public :: simtype
   PRIVATE
   SEQUENCE
   !! Dielectric Constant for reaction field
   real ( kind = dp ) :: dielect = 0.0_dp
   !! Periodic Boundry Conditions
   logical :: SIM_uses_PBC
   logical :: SIM_uses_LJ
   logical :: SIM_uses_sticky
   logical :: SIM_uses_charges
   logical :: SIM_uses_dipoles
   logical :: SIM_uses_RF
   logical :: SIM_uses_GB
   logical :: SIM_uses_EAM
end type simtype
#endif
