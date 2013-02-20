/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
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
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *                                                                      
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "config.h"
#ifdef  IS_MPI
#include <mpi.h>
#endif

int nChecks;

#include "utils/simError.h"

errorStruct painCave;

char checkPointMsg[MAX_SIM_ERROR_MSG_LENGTH];
int worldRank;

void initSimError( void ){
  painCave.errMsg[0] = '\0';
  painCave.isFatal = 0;
  painCave.severity = OPENMD_ERROR;
  painCave.isEventLoop = 0;
  nChecks = 0;
#ifdef IS_MPI
  worldRank = MPI::COMM_WORLD.Get_rank();
#else
  worldRank = 0;
#endif
}

int simError( void ) {
  
  char errorMsg[MAX_SIM_ERROR_MSG_LENGTH];

#ifdef IS_MPI
  int myError = 1;
  int isError;
  char nodeMsg[MAX_SIM_ERROR_MSG_LENGTH];
#endif
  
  strcpy(errorMsg, "OpenMD ");
  switch( painCave.severity ) {
  case OPENMD_WARNING:
    strcat(errorMsg, "warning");
    break;
  case OPENMD_INFO:
    strcat(errorMsg, "info");
    break;
  default:
    if( painCave.isFatal ) {
      strcat(errorMsg, "FATAL ");
    }
    strcat(errorMsg, "ERROR");
  }
  
#ifdef IS_MPI
  if (worldRank == 0) {
    if ( painCave.isEventLoop ) {
      sprintf( nodeMsg, " (reported by MPI node %d)", worldRank);
      strncat(errorMsg, nodeMsg, strlen(nodeMsg));
    }
#endif
    
    strcat(errorMsg, ":\n\t");
    
    strncat(errorMsg, painCave.errMsg, strlen(painCave.errMsg));
    
    strcat(errorMsg, "\n");
    fprintf(stderr, "%s", errorMsg);
    
#ifdef IS_MPI
    if (painCave.isEventLoop) 
      return 1;
  }
#endif   

  if (painCave.isFatal) {
#ifdef IS_MPI    
    MPI::COMM_WORLD.Allreduce(&myError, &isError, 1, MPI::INT, MPI::LOR);
    MPI::Finalize();
#endif
    exit(0);
  }  
  return 1;  
}


void errorCheckPoint( void ){
  
  int myError = 0;
  int isError = 0;
  
#ifdef IS_MPI
  MPI::COMM_WORLD.Allreduce(&myError, &isError, 1, MPI::INT, MPI::LOR);
#else
  isError = myError;
#endif
  
  if( isError ){    
#ifdef IS_MPI
    MPI::Finalize();
#endif    
    exit(0);
  }
  
#ifdef CHECKPOINT_VERBOSE  
  nChecks++;

#ifdef IS_MPI
  if( worldRank == 0 ){
#endif

    fprintf(stderr,
	    "Checkpoint #%d reached: %s\n",
	    nChecks,
	    checkPointMsg );
#ifdef IS_MPI
  }
#endif

#endif 
}


