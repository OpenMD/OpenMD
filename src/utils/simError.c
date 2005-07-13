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
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "config.h"

#ifdef  IS_MPI
#include <mpi.h>

int nChecks;
#endif 

#include "utils/simError.h"

errorStruct painCave;

#ifdef IS_MPI

char checkPointMsg[MAX_SIM_ERROR_MSG_LENGTH];
int worldRank;

#endif


void initSimError( void ){
  painCave.errMsg[0] = '\0';
  painCave.isFatal = 0;
  painCave.severity = OOPSE_ERROR;
#ifdef IS_MPI
  painCave.isEventLoop = 0;
  nChecks = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &worldRank );
#endif
}

int simError( void ) {
  
  int myError = 1;
  int isError;
  char errorMsg[MAX_SIM_ERROR_MSG_LENGTH];
  char nodeMsg[MAX_SIM_ERROR_MSG_LENGTH];
  
  strcpy(errorMsg, "OOPSE ");
  switch( painCave.severity ) {
  case OOPSE_WARNING:
    strcat(errorMsg, "WARNING");
    break;
  case OOPSE_INFO:
    strcat(errorMsg, "INFO");
    break;
  default:
    if( painCave.isFatal ) {
      strcat(errorMsg, "FATAL ");
    }
    strcat(errorMsg, "ERROR");
  }
  
#ifdef IS_MPI
  if ( painCave.isEventLoop ) {
    sprintf( nodeMsg, " (reported by MPI node %d)", worldRank);
    strncat(errorMsg, nodeMsg, strlen(nodeMsg));
  }
#endif

  strcat(errorMsg, ":\n\t");

  strncat(errorMsg, painCave.errMsg, strlen(painCave.errMsg));

  strcat(errorMsg, "\n");
  fprintf(stderr, errorMsg);

#ifdef IS_MPI
  if (painCave.isEventLoop) 
    return 1;
#endif   

  if (painCave.isFatal) {
#ifdef IS_MPI    
    MPI_Allreduce( &myError, &isError, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD );
    MPI_Finalize();
#endif
    exit(0);
  } 

  return 1;  
}
 
  
#ifdef IS_MPI

void MPIcheckPoint( void ){
  
  int myError = 0;
  int isError;

  MPI_Allreduce( &myError, &isError, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD );
  if( isError ){
    MPI_Finalize();
    exit(0);
  }

#ifdef CHECKPOINT_VERBOSE  
  nChecks++;
  if( worldRank == 0 ){
    fprintf(stderr,
	    "Checkpoint #%d reached: %s\n",
	    nChecks,
	    checkPointMsg );
  }
#endif 

}

#endif 
