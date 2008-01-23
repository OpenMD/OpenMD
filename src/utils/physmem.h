#ifndef UTILS_PHYSMEM_H
#define UTILS_PHYSMEM_H

#ifdef  __cplusplus
extern "C" {
#endif

  /** Return the total amount of physical memory.  */
  double physmem_total (void);

  /** Return the amount of physical memory available.  */
  double physmem_available (void);

#ifdef  __cplusplus
}
#endif

#endif
