#ifdef __C
#ifndef __FCOULOMBICCORRECTION
#define __FCOULOMBICCORRECTION

#define NONE                   1
#define UNDAMPED_WOLF          2
#define WOLF                   3
#define REACTION_FIELD         4

#endif
#endif /*__C*/

#ifdef  __FORTRAN90

  INTEGER, PARAMETER:: NONE                      = 1
  INTEGER, PARAMETER:: UNDAMPED_WOLF             = 2
  INTEGER, PARAMETER:: WOLF                      = 3
  INTEGER, PARAMETER:: REACTION_FIELD            = 4

#endif
