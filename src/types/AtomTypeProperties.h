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
  int is_Electrostatic;
  int is_Charge;
  int is_Dipole;
  int is_Sticky;
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
   logical :: is_Directional
   logical :: is_LennardJones
   logical :: is_Electrostatic
   logical :: is_Charge
   logical :: is_Dipole
   logical :: is_Sticky
   logical :: is_GayBerne
   logical :: is_EAM
   logical :: is_Shape
   logical :: is_FLARB
end type AtomTypeProperties
#endif
