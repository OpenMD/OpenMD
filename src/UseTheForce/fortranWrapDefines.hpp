#ifndef __FORTRAN_WRAP_DEFINES_H__
#define __FORTRAN_WRAP_DEFINES_H__

#define __C
#include "brains/fSimulation.h"

// here we declare the function pointer typedefs for fortran functions

extern "C" {

  typedef void (*makeAtype_TD) ( int* unique_ident, 
				 int* isLJ, 
				 int* isSticky,
				 int* isDipole, 
				 int* isGB, 
				 int* isEAM,
                                 int* isCharge,
				 double* lj_epslon, 
				 double* lj_sigma, 
                                 double* charge,
				 double* dipole_moment, 
				 int* status );

  typedef void (*newEAMtype_TD)( double* lattice_constant, 
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
  
  typedef void (*initFortranFF_TD)( int* LJ_mix_policy, 
				    int* useReactionField,
				    int *isError );
  
  typedef void (*doForceLoop_TD)( double* positionArray,
  				  double* rcArray,
				  double* RotationMatrixArray,
				  double* unitVectorArray_l,
				  double* forceArray,
				  double *torqueArray,
				  double* StressTensor, 
				  double* potentialEnergy, 
				  short int* doPotentialCalc, 
				  short int* doStressCalc,
				  int* isError );
				 
  typedef void (*set_sticky_params_TD)( double* sticky_w0, 
					double* sticky_v0,
					double* sticky_v0p,
					double* sticky_rl,
					double* sticky_ru,
					double* sticky_rlp,
					double* sticky_rup );

  typedef void (*set_gb_pair_params_TD)( double* GB_sigma,
					 double* GB_l2b_ratio,
					 double* GB_eps,
					 double* GB_eps_ratio,
					 double* GB_mu,
					 double* GB_nu );

  typedef void (*setFortranSim_TD)( simtype* the_Info,
				    int* nGlobal, 
				    int* nLocal, 
				    int* identArray,
				    int* nLocalExcludes,
				    int* excludesLocalArray,
				    int* nGlobalExcludes,
				    int* excludesGlobalArray,
				    int* molMembershipArray,
				    double* mfact,
				    int* ngroup,
				    int* globalGroupMembership,
				    int* isError );

  typedef void (*setFortranBox_TD) ( double *Hmat,
				     double *HmatI,
				     int* orthoRhombic );

  typedef void (*notifyFortranCutOff_TD) ( double *rCut,
                                           double *rSw,
					   double *rList );
}


#ifdef IS_MPI

#include "UseTheForce/mpiComponentPlan.h"



extern "C" {
  
  typedef void (*setFortranMPI_TD)( mpiSimData* the_mpiPlug,
				    int* nLocal, 
				    int* globalAtomIndex,
                                    int* nGroupsLocal,
                                    int* globalGroupIndex,
				    int* isError );

}

#endif // is_mpi

#endif // frotranWrapDefines.hpp
