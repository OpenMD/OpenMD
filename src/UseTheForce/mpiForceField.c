#ifdef IS_MPI

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "mpiForceField.h"




void sendFrcStruct( void *frcStruct, MPI_Datatype structType ){
  
  MPI_Bcast(frcStruct,1,structType,0,MPI_COMM_WORLD);
}


void receiveFrcStruct( void *frcStruct, MPI_Datatype structType ){
  
  MPI_Bcast(frcStruct,1,structType,0,MPI_COMM_WORLD);
}


#endif // is_mpi
