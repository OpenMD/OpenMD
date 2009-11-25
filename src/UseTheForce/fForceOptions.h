#ifdef __OPENMD_C
#ifndef __FFORCEOPTIONS
#define __FFORCEOPTIONS

#define GEOMETRIC_MIXING_RULE  1
#define ARITHMETIC_MIXING_RULE 2
#define CUBIC_MIXING_RULE 3
#define HHG_MIXING_RULE 4

typedef  struct{
  int DistanceMixingRule;
  int EnergyMixingRule;
  RealType vdw12scale;
  RealType vdw13scale;
  RealType vdw14scale;
  RealType electrostatic12scale;  
  RealType electrostatic13scale;  
  RealType electrostatic14scale;  
  RealType GayBerneMu;  
  RealType GayBerneNu;  
} ForceOptions;


#endif
#endif /*__OPENMD_C*/

#ifdef  __FORTRAN90

  INTEGER, PARAMETER:: GEOMETRIC_MIXING_RULE  = 1
  INTEGER, PARAMETER:: ARITHMETIC_MIXING_RULE = 2
  INTEGER, PARAMETER:: CUBIC_MIXING_RULE = 3
  INTEGER, PARAMETER:: HHG_MIXING_RULE = 4

  type :: ForceOptions
    SEQUENCE
    integer :: DistanceMixingRule
    integer :: EnergyMixingRule
    real(kind=dp) :: vdw12scale
    real(kind=dp) :: vdw13scale
    real(kind=dp) :: vdw14scale
    real(kind=dp) :: electrostatic12scale
    real(kind=dp) :: electrostatic13scale
    real(kind=dp) :: electrostatic14scale
    real(kind=dp) :: GayBerneMu
    real(kind=dp) :: GayBerneNu
  end type ForceOptions

#endif
