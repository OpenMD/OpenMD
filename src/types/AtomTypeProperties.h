 /*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 */
 
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
  int is_Quadrupole;
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
   integer :: is_Directional
   integer :: is_LennardJones
   integer :: is_Charge
   integer :: is_Dipole
   integer :: is_Quadrupole
   integer :: is_Sticky
   integer :: is_GayBerne
   integer :: is_EAM
   integer :: is_Shape
   integer :: is_FLARB
end type AtomTypeProperties
#endif
