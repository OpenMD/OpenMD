#ifndef UTILS_PHYSMEM_H
#define UTILS_PHYSMEM_H

#ifdef  __cplusplus
extern "C" {
#endif

  /** Return the total amount of physical memory.  */
  RealType physmem_total ();

  /** Return the amount of physical memory available.  */
  RealType physmem_available ();

#ifdef  __cplusplus
}
#endif

#endif
