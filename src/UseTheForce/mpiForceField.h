#ifndef __MPIFORCEFIELD_H__
#define __MPIFORCEFIELD_H__

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

  void sendFrcStruct( void *frcStruct, MPI_Datatype structType );

  void receiveFrcStruct( void *frcStruct, MPI_Datatype structType );

#ifdef __cplusplus
}
#endif

#endif
