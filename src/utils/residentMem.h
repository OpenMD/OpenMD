#ifndef UTILS_RESIDENTMEM_H
#define UTILS_RESIDENTMEM_H

#ifdef  __cplusplus
extern "C" {
#endif

/** Returns an estimate of the amount of memory being used by other processes.  */
double residentMem ();

#ifdef  __cplusplus
}
#endif

#endif
