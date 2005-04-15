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
 
#ifndef TYPES_ATOMSTAMP_HPP
#define TYPES_ATOMSTAMP_HPP

#include "io/LinkedAssign.hpp"

class AtomStamp{
  
 public:
  AtomStamp();
  ~AtomStamp();

  void setPosition( double x, double y, double z );
  void setOrientation( double phi, double theta, double psi );
  char* assignString( char* lhs, char* rhs );
  char* assignDouble( char* lhs, double rhs );
  char* assignInt( char* lhs, int rhs );
  char* checkMe( void );

  char* getType( void ) { return type; }
  short int havePosition( void ) { return have_position; }
  short int haveOrientation( void ) { return have_orientation; }
  double getPosX( void ) { return pos[0]; }
  double getPosY( void ) { return pos[1]; }
  double getPosZ( void ) { return pos[2]; }
  double getEulerPhi( void )   { return ornt[0]; }
  double getEulerTheta( void ) { return ornt[1]; }
  double getEulerPsi( void )   { return ornt[2]; }
  

 private:

  double pos[3]; //the position vector
  short int have_position; // boolean for positions
  double ornt[3]; // the Euler angles
  short int have_orientation;
  char type[100]; // the type name of the atom
  short int have_type;
  
  LinkedAssign* unhandled; // the list of unhandled assignments
  short int have_extras;
};

#endif
