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
 
#ifndef __FORTRAN90
#ifndef UTILS_SIMERROR_H
#define UTILS_SIMERROR_H

#define MAX_SIM_ERROR_MSG_LENGTH 2000

#define OOPSE_ERROR   1
#define OOPSE_WARNING 2
#define OOPSE_INFO    3

typedef struct{
  char errMsg[MAX_SIM_ERROR_MSG_LENGTH];
  int isFatal;
  int severity;
#ifdef IS_MPI
  int isEventLoop;
#endif // IS_MPI
} errorStruct;

extern errorStruct painCave;

#ifdef IS_MPI

extern char checkPointMsg[MAX_SIM_ERROR_MSG_LENGTH];

extern int worldRank;
#endif

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus
  
  int simError( void ); // returns 1 if handled. 0 otherwise.

  void initSimError( void ); // needed to be called from main before anything
                             // goes wrong.

#ifdef IS_MPI
  
  void MPIcheckPoint( void );
  
#endif // IS_MPI

#ifdef __cplusplus
}
#endif //__cplusplus

#endif // __SIMERROR_H__

#else // __FORTRAN90

  INTEGER, PARAMETER:: OOPSE_ERROR   = 1
  INTEGER, PARAMETER:: OOPSE_WARNING = 2
  INTEGER, PARAMETER:: OOPSE_INFO    = 3
  INTEGER, PARAMETER:: MAX_SIM_ERROR_MSG_LENGTH = 2000
  
type, public :: errorStruct
  PRIVATE
  SEQUENCE
  character(len = MAX_SIM_ERROR_MSG_LENGTH) :: errMsg
  logical :: isFatal
  integer :: severity
#ifdef IS_MPI
  logical :: isEventLoop;
#endif // IS_MPI
end type errorStruct

type (errorStruct), public, save :: painCave

#endif // __FORTRAN90
