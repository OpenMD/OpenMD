#ifdef __OPENMD_C

#ifndef __FSIMULATION

#define __FSIMULATION
/** This header provides dual access for the simulation structure between 
    fortran and C for the simtype structure. NOTE: Sequence of struct 
    components must match between C and fortran and in general be packed 
    RealType,int,char. 
*/
typedef  struct{
  int SIM_uses_PBC;
  int SIM_uses_DirectionalAtoms;
  int SIM_uses_MetallicAtoms;
  int SIM_uses_AtomicVirial;
  int SIM_requires_SkipCorrection;
  int SIM_requires_SelfCorrection;
} simtype;
#endif /*__FSIMULATION*/
#endif /*__OPENMD_C*/

#ifdef  __FORTRAN90

  type, public :: simtype
    PRIVATE
    SEQUENCE
    !! Periodic Boundry Conditions
    logical :: SIM_uses_PBC
    logical :: SIM_uses_DirectionalAtoms
    logical :: SIM_uses_MetallicAtoms
    logical :: SIM_uses_AtomicVirial
    logical :: SIM_requires_SkipCorrection
    logical :: SIM_requires_SelfCorrection
  end type simtype

#endif
  
