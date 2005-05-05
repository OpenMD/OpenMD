#ifdef __C
#ifndef TYPES_ATOMTYPEPROPERTIES_H
#define TYPES_ATOMTYPEPROPERTIES_H

/** 
 * This header provides dual access for the AtomTypeProperties between
 * fortran and C. NOTE: The sequence of struct components MUST match
 * between C and Fortran and in general be packed double,int,char.
 */
typedef  struct{
  int ident;
  int is_Directional;
  int is_LennardJones;
  int is_Charge;
  int is_Dipole;
  int is_SplitDipole;
  int is_Quadrupole;
  int is_Sticky;
  int is_StickyPower;
  int is_GayBerne;
  int is_EAM;
  int is_Shape;
  int is_FLARB;
} AtomTypeProperties;
#endif 
#endif 

#ifdef  __FORTRAN90

  type :: AtomTypeProperties
    SEQUENCE
    integer :: ident
    integer :: is_Directional
    integer :: is_LennardJones
    integer :: is_Charge
    integer :: is_Dipole
    integer :: is_SplitDipole
    integer :: is_Quadrupole
    integer :: is_Sticky
    integer :: is_StickyPower
    integer :: is_GayBerne
    integer :: is_EAM
    integer :: is_Shape
    integer :: is_FLARB
  end type AtomTypeProperties
#endif
