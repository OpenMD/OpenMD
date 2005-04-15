#ifndef UTILS_PHYSMEM_H
#define UTILS_PHYSMEM_H

#ifdef  __cplusplus
extern "C" {
#endif

  /** Return the total amount of physical memory.  */
  double physmem_total ();

  /** Return the amount of physical memory available.  */
  double physmem_available ();

#ifdef  __cplusplus
}
#endif

#endif
