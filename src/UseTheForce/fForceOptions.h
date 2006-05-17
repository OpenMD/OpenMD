#ifdef __C
#ifndef __FFORCEOPTIONS
#define __FFORCEOPTIONS

#define GEOMETRIC_MIXING_RULE  1
#define ARITHMETIC_MIXING_RULE 2

typedef  struct{
  int DistanceMixingRule;
  int EnergyMixingRule;
  RealType vdw14scale;
  RealType electrostatic14scale;  
} ForceOptions;


#endif
#endif /*__C*/

#ifdef  __FORTRAN90

  INTEGER, PARAMETER:: GEOMETRIC_MIXING_RULE  = 1
  INTEGER, PARAMETER:: ARITHMETIC_MIXING_RULE = 2

  type :: ForceOptions
    SEQUENCE
    integer :: DistanceMixingRule
    integer :: EnergyMixingRule
    real(kind=dp) :: vdw14scale
    real(kind=dp) :: electrostatic14scale
  end type ForceOptions

#endif
