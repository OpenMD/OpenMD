#ifndef __MPIBASS_H
#define __MPIBASS_H
#include "config.h"

#define MPI_INTERFACE_ABORT 2
#define MPI_INTERFACE_DONE  1
#define MPI_INTERFACE_CONTINUE 0

#include <mpi.h>

#ifndef __is_lex__
#include "BASS_interface.h"
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

// Structure to pass mpi a BASS event

typedef struct mpiBASSEvent{
  int type;
  double d1,d2,d3;
  int    i1;
  char   cArray[120];
  char   lhs[80];
} mBEvent;

// types for mpiBASSEvent.type
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


// Define the mpi datatype
#ifdef __mpiBASSEVENT
MPI_Datatype mpiBASSEventType;
#endif

#endif
