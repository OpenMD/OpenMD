
#ifndef EAM_MODULE_INTERFACE_H
#define EAM_MODULE_INTERFACE_H

#define __C
#include "config.h"
extern "C"{
  void F90_FUNC(neweamtype, NEWEAMTYPE)( double* lattice_constant, 
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
                   int* status ){
    
    F90_FUNC(neweamtype, NEWEAMTYPE)(lattice_constant, 
                                     eam_nrho,
                                     eam_drho,
                                     eam_nr,
                                     eam_dr,
                                     eam_rcut,
                                     eam_rvals, 
                                     eam_rhovals,
                                     eam_Frhovals,
                                     eam_ident,
                                     status );
  }
}
#endif