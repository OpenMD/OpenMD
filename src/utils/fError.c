/* Provides a fortran - c interface for info writer system.
 */
#include <string.h>
#include "config.h"

#include "simError.h" 

void F90_FUNC_(c_simerror, C_SIMERROR) (errorStruct* pc);

void F90_FUNC_(c_simerror, C_SIMERROR) (errorStruct* pc){

  strcpy(painCave.errMsg, (*pc).errMsg);
  painCave.severity = (*pc).severity;
  painCave.isFatal = (*pc).isFatal;
#ifdef IS_MPI
  painCave.isEventLoop = (*pc).isEventLoop;
#endif
  simError();

}
