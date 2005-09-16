#ifdef __C
#ifndef __FELECTROSTATICSUMMATIONMETHOD
#define __FELECTROSTATICSUMMATIONMETHOD

#define NONE                   1
#define UNDAMPED_WOLF          2
#define DAMPED_WOLF            3
#define REACTION_FIELD         4
/* Ewald methods aren't supported yet */
#define EWALD_FULL             5 
#define EWALD_PME              6
#define EWALD_SPME             7

#endif
#endif /*__C*/

#ifdef  __FORTRAN90

  INTEGER, PARAMETER:: NONE                      = 1
  INTEGER, PARAMETER:: UNDAMPED_WOLF             = 2
  INTEGER, PARAMETER:: DAMPED_WOLF               = 3
  INTEGER, PARAMETER:: REACTION_FIELD            = 4
  INTEGER, PARAMETER:: EWALD_FULL                = 5
  INTEGER, PARAMETER:: EWALD_PME                 = 6
  INTEGER, PARAMETER:: EWALD_SPME                = 7

#endif
