#ifdef __OPENMD_C
#ifndef __FSWITCHINGFUNCTIONTYPE
#define __FSWITCHINGFUNCTIONTYPE

#define CUBIC                1
#define FIFTH_ORDER_POLY     2

#endif
#endif /*__OPENMD_C*/

#ifdef  __FORTRAN90

  INTEGER, PARAMETER:: CUBIC                = 1
  INTEGER, PARAMETER:: FIFTH_ORDER_POLY     = 2

#endif
