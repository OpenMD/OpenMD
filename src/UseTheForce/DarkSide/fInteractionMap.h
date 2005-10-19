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

#define LR_POT_TYPES      4

#define VDW_POT           0
#define ELECTROSTATIC_POT 1
#define HB_POT            2
#define METALLIC_POT      3


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

  INTEGER, PARAMETER:: LR_POT_TYPES      =   4
  !! remember fortran is shifted up by one for arrays
  INTEGER, PARAMETER:: VDW_POT           =   1
  INTEGER, PARAMETER:: ELECTROSTATIC_POT =   2
  INTEGER, PARAMETER:: HB_POT            =   3
  INTEGER, PARAMETER:: METALLIC_POT      =   4

#endif
