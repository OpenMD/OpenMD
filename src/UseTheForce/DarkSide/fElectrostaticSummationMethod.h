#ifdef __C
#ifndef __FELECTROSTATICSUMMATIONMETHOD
#define __FELECTROSTATICSUMMATIONMETHOD

#define NONE                   1
#define SWITCHING_FUNCTION     2
#define SHIFTED_POTENTIAL      3
#define SHIFTED_FORCE          4
#define REACTION_FIELD         5
/* Ewald methods aren't supported yet */
#define EWALD_FULL             6 
#define EWALD_PME              7
#define EWALD_SPME             8

#endif
#endif /*__C*/

#ifdef  __FORTRAN90

  INTEGER, PARAMETER:: NONE                      = 1
  INTEGER, PARAMETER:: SWITCHING_FUNCTION        = 2
  INTEGER, PARAMETER:: SHIFTED_POTENTIAL         = 3
  INTEGER, PARAMETER:: SHIFTED_FORCE             = 4
  INTEGER, PARAMETER:: REACTION_FIELD            = 5
  INTEGER, PARAMETER:: EWALD_FULL                = 6
  INTEGER, PARAMETER:: EWALD_PME                 = 7
  INTEGER, PARAMETER:: EWALD_SPME                = 8

#endif
