#ifndef TYPES_ATOMTYPEPROPERTIES_H
#define TYPES_ATOMTYPEPROPERTIES_H

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
  int is_SC;
} AtomTypeProperties;
#endif 
