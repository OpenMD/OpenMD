#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "config.h"

#ifdef  IS_MPI
#include <mpi.h>

int nChecks;
#endif // IS_MPI

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
#endif // CHECKPOINT_VERBOSE  

}

#endif // IS_MPI
