#ifndef __FORTRAN90
#ifndef __SIMERROR_H__
#define __SIMERROR_H__

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
