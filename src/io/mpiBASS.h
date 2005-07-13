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
 
#ifndef IO_MPIBASS_H
#define IO_MPIBASS_H
#include "config.h"

#define MPI_INTERFACE_ABORT 2
#define MPI_INTERFACE_DONE  1
#define MPI_INTERFACE_CONTINUE 0

#include <mpi.h>

#ifndef __is_lex__
#include "io/BASS_interface.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

  void mpiInterfaceExit(void);
  void mpiEventInit(void);
#ifndef __is_lex__
  void throwMPIEvent(event* event);
  void mpiEventLoop(void);
#endif

#ifdef __cplusplus
}
#endif

/* Structure to pass mpi a BASS event*/

typedef struct mpiBASSEvent{
  int type;
  double d1,d2,d3;
  int    i1;
  char   cArray[120];
  char   lhs[80];
} mBEvent;

/* types for mpiBASSEvent.type*/
#define mpiMOLECULE     0
#define mpiATOM         1
#define mpiBOND         2
#define mpiBEND         3
#define mpiTORSION      4
#define mpiCOMPONENT    5
#define mpiPOSITION     6
#define mpiASSIGNMENT_i 7
#define mpiASSIGNMENT_d 8
#define mpiASSIGNMENT_s 9
#define mpiMEMBERS      10
#define mpiCONSTRAINT   11
#define mpiORIENTATION  12
#define mpiBLOCK_END    13
#define mpiZCONSTRAINT  14
#define mpiRIGIDBODY    15
#define mpiCUTOFFGROUP  16
#define mpiMEMBER       17


/* Define the mpi datatype*/
#ifdef __mpiBASSEVENT
MPI_Datatype mpiBASSEventType;
#endif

#endif
