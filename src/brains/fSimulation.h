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

#ifndef __FSIMULATION

#define __FSIMULATION
/** This header provides dual access for the simulation structure between 
    fortran and C for the simtype structure. NOTE: Sequence of struct 
    components must match between C and fortran and in general be packed 
    double,int,char. 
*/
typedef  struct{
  double dielect;
  int SIM_uses_PBC;
  int SIM_uses_DirectionalAtoms;
  int SIM_uses_LennardJones;
  int SIM_uses_Electrostatics;
  int SIM_uses_Charges;
  int SIM_uses_Dipoles;
  int SIM_uses_Sticky;
  int SIM_uses_GayBerne;
  int SIM_uses_EAM;
  int SIM_uses_Shapes;
  int SIM_uses_FLARB;
  int SIM_uses_RF;
} simtype;
#endif //__FSIMULATION
#endif //__C

#ifdef  __FORTRAN90

type, public :: simtype
   PRIVATE
   SEQUENCE
   !! Dielectric Constant for reaction field
   real ( kind = dp ) :: dielect = 0.0_dp
   !! Periodic Boundry Conditions
   logical :: SIM_uses_PBC
   logical :: SIM_uses_DirectionalAtoms
   logical :: SIM_uses_LennardJones
   logical :: SIM_uses_Electrostatics
   logical :: SIM_uses_Charges
   logical :: SIM_uses_Dipoles
   logical :: SIM_uses_Sticky
   logical :: SIM_uses_GayBerne
   logical :: SIM_uses_EAM
   logical :: SIM_uses_Shapes
   logical :: SIM_uses_FLARB
   logical :: SIM_uses_RF
end type simtype
#endif

