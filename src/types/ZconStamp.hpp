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
 
#ifndef TYPES_ZCONSTAMP_HPP
#define TYPES_ZCONSTAMP_HPP

#include "io/LinkedAssign.hpp"

class ZconStamp{

public:
  ZconStamp(int theIndex);
  ~ZconStamp();

  int assignString( char* lhs, char* rhs, char** err );
  int assignDouble( char* lhs, double rhs, char** err );
  int assignInt( char* lhs, int rhs, char** err );

  int haveExtras( void ) { return have_extras; }
  LinkedAssign* getExtras( void ) { return unhandled; }

  char* checkMe( void );

  int    getMolIndex( void ) { return molIndex; }
  double getZpos( void )     { return zPos; }
  double getKratio( void )   { return kRatio; }
  double getCantVel( void )   { return cantVel; }
  
  short int haveZpos( void ) { return have_zPos; }
  short int haveKratio( void ) { return have_kRatio; }
  short int haveCantVel( void ) { return have_cantVel; }
  
private:

  short int have_zPos;
  short int have_molIndex;
  short int have_kRatio;
  short int have_cantVel;

  double zPos;
  double kRatio;
  double cantVel;
  int molIndex;
  
  int index;

  LinkedAssign* unhandled; // the unhandled assignments
  short int have_extras;  
};


#endif // __ZCONSTAMP_H__
