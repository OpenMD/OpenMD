#ifndef UTILS_PHYSMEM_H
#define UTILS_PHYSMEM_H

extern "C" {
/** Return the total amount of physical memory.  */
double physmem_total ();

/** Return the amount of physical memory available.  */
double physmem_available ();
}
#endif
