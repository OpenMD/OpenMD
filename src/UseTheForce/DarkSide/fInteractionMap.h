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

#endif
