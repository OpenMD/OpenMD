#ifndef USETHEFORCE_DARKSIDE_EAM_INTERFACE_H
#define USETHEFORCE_DARKSIDE_EAM_INTERFACE_H

#define __C

#include "config.h"

#define newEAMtype F90_FUNC(neweamtype, NEWEAMTYPE)

extern "C"{
  void newEAMtype( double* lattice_constant, 
                   int* eam_nrho,
                   double* eam_drho,
                   int* eam_nr,
                   double* eam_dr,
                   double* eam_rcut,
                   double* eam_rvals, 
                   double* eam_rhovals,
                   double* eam_Frhovals,
                   int* eam_ident,
                   int* status );
}  
#endif
