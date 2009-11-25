 
#ifndef __FORTRAN90
#ifndef UTILS_SIMERROR_H
#define UTILS_SIMERROR_H

#define MAX_SIM_ERROR_MSG_LENGTH 2000

#define OPENMD_ERROR   1
#define OPENMD_WARNING 2
#define OPENMD_INFO    3

typedef struct{
  char errMsg[MAX_SIM_ERROR_MSG_LENGTH];
  int isFatal;
  int severity;
  int isEventLoop;
} errorStruct;

extern errorStruct painCave;

extern char checkPointMsg[MAX_SIM_ERROR_MSG_LENGTH];

extern int worldRank;

#ifdef __cplusplus
extern "C" {
#endif 
  
  int simError( void ); 

  void initSimError( void ); 

  void errorCheckPoint( void );
                             
#ifdef __cplusplus
}
#endif 

#endif 

#else 

  INTEGER, PARAMETER:: OPENMD_ERROR   = 1
  INTEGER, PARAMETER:: OPENMD_WARNING = 2
  INTEGER, PARAMETER:: OPENMD_INFO    = 3
  INTEGER, PARAMETER:: MAX_SIM_ERROR_MSG_LENGTH = 2000
  
  type, public :: errorStruct
    PRIVATE
    SEQUENCE
    character(len = MAX_SIM_ERROR_MSG_LENGTH) :: errMsg
    logical :: isFatal
    integer :: severity
    logical :: isEventLoop;
  end type errorStruct

  type (errorStruct), public, save :: painCave

#endif 
