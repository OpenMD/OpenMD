#ifdef __C

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
  int SIM_uses_LennardJones;
  int SIM_uses_Electrostatics;
  int SIM_uses_Charges;
  int SIM_uses_Dipoles;
  int SIM_uses_Sticky;
  int SIM_uses_StickyPower;
  int SIM_uses_GayBerne;
  int SIM_uses_EAM;
  int SIM_uses_Shapes;
  int SIM_uses_FLARB;
  int SIM_uses_RF;
  int SIM_uses_SF;
  int SIM_uses_SP;
  int SIM_uses_SC;
  int SIM_uses_MEAM;
  int SIM_uses_BoxDipole;
  int SIM_uses_AtomicVirial;
} simtype;
#endif /*__FSIMULATION*/
#endif /*__C*/

#ifdef  __FORTRAN90

  type, public :: simtype
    PRIVATE
    SEQUENCE
    !! Periodic Boundry Conditions
    logical :: SIM_uses_PBC
    logical :: SIM_uses_DirectionalAtoms
    logical :: SIM_uses_LennardJones
    logical :: SIM_uses_Electrostatics
    logical :: SIM_uses_Charges
    logical :: SIM_uses_Dipoles
    logical :: SIM_uses_Sticky
    logical :: SIM_uses_StickyPower
    logical :: SIM_uses_GayBerne
    logical :: SIM_uses_EAM
    logical :: SIM_uses_Shapes
    logical :: SIM_uses_FLARB
    logical :: SIM_uses_RF
    logical :: SIM_uses_SF
    logical :: SIM_uses_SP
    logical :: SIM_uses_SC
    logical :: SIM_uses_MEAM
    logical :: SIM_uses_BoxDipole
    logical :: SIM_uses_AtomicVirial
  end type simtype

#endif
  
