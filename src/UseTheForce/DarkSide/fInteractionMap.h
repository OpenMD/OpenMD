#ifdef __C
#ifndef __FINTERACTIONMAP
#define __FINTERACTIONMAP

#define LJ_PAIR              1
#define ELECTROSTATIC_PAIR   2
#define STICKY_PAIR          4
#define STICKYPOWER_PAIR     8
#define GB_PAIR             16
#define GB_LJ               32
#define EAM_PAIR            64
#define SHAPE_PAIR         128
#define SHAPE_LJ           256
#define FLARB_PAIR         512

#define POT_ARRAY_SIZE 11

#define LJ_POT            1
#define ELECTROSTATIC_POT 2
#define STICKY_POT        3
#define STICKYPOWER_POT   4
#define GB_POT            5
#define GB_LJ_POT         6
#define EAM_POT           7
#define SHAPE_POT         8
#define SHAPE_LJ_POT      9
#define FLARB_POT         10
#define RF_POT            11

#endif
#endif /*__C*/

#ifdef  __FORTRAN90

  INTEGER, PARAMETER:: LJ_PAIR            =   1
  INTEGER, PARAMETER:: ELECTROSTATIC_PAIR =   2
  INTEGER, PARAMETER:: STICKY_PAIR        =   4
  INTEGER, PARAMETER:: STICKYPOWER_PAIR   =   8
  INTEGER, PARAMETER:: GAYBERNE_PAIR      =  16
  INTEGER, PARAMETER:: GAYBERNE_LJ        =  32
  INTEGER, PARAMETER:: EAM_PAIR           =  64
  INTEGER, PARAMETER:: SHAPE_PAIR         = 128
  INTEGER, PARAMETER:: SHAPE_LJ           = 256
  INTEGER, PARAMETER:: FLARB_PAIR         = 512

  INTEGER, PARAMETER:: POT_ARRAY_SIZE    =   10

  INTEGER, PARAMETER:: LJ_POT            =   1
  INTEGER, PARAMETER:: ELECTROSTATIC_POT =   2
  INTEGER, PARAMETER:: STICKY_POT        =   3
  INTEGER, PARAMETER:: STICKYPOWER_POT   =   4
  INTEGER, PARAMETER:: GAYBERNE_POT      =   5
  INTEGER, PARAMETER:: GAYBERNE_LJ_POT   =   6
  INTEGER, PARAMETER:: EAM_POT           =   7
  INTEGER, PARAMETER:: SHAPE_POT         =   8
  INTEGER, PARAMETER:: SHAPE_LJ_POT      =   9
  INTEGER, PARAMETER:: FLARB_POT         =   10
  INTEGER, PARAMETER:: RF_POT            =   11



#endif
