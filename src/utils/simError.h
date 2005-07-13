 
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
#endif 
} errorStruct;

extern errorStruct painCave;

#ifdef IS_MPI

extern char checkPointMsg[MAX_SIM_ERROR_MSG_LENGTH];

extern int worldRank;
#endif

#ifdef __cplusplus
extern "C" {
#endif 
  
  int simError( void ); 

  void initSimError( void ); 
                             

#ifdef IS_MPI
  
  void MPIcheckPoint( void );
  
#endif 

#ifdef __cplusplus
}
#endif 

#endif 

#else 

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
#endif 
  end type errorStruct

  type (errorStruct), public, save :: painCave

#endif 
