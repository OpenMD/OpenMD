#define __C

#include "config.h"
#include "brains/fSimulation.h"
#include "UseTheForce/fortranWrappers.hpp"

// declare the actual instances of the function pointers

makeAtype_TD makeAtype;
initFortranFF_TD initFortranFF;
set_sticky_params_TD set_sticky_params;
set_gb_pair_params_TD set_gb_pair_params;
newEAMtype_TD newEAMtype;

// declare the functions on the fortran side

extern "C" {
  
  typedef void (*ffWrapFunction_TD)( makeAtype_TD p1,
				     initFortranFF_TD p2,
				     doForceLoop_TD p3,
				     set_sticky_params_TD p4,
				     set_gb_pair_params_TD p5,
				     newEAMtype_TD p6 );

  void F90_FUNC(wrapforcefield, WRAPFORCEFIELD)( ffWrapFunction_TD myWF );
  
  typedef void (*smWrapFunction_TD)( setFortranSim_TD p1,
				     setFortranBox_TD p2,
				     notifyFortranCutOff_TD p3 );

  void F90_FUNC(wrapsimmod, WRAPSIMMOD) ( smWrapFunction_TD myWF );
  
#ifdef IS_MPI
  
  typedef void (*spmWrapFunction_TD)( setFortranMPI_TD p1 );

  void F90_FUNC(wrapsimparallelmod, WRAPSIMPARALLELMOD)(spmWrapFunction_TD myWF);

#endif // is_mpi
}

// declare the functions that are defined in this file
extern "C"{
  void wrapFF(makeAtype_TD p1,
	      initFortranFF_TD p2,
	      doForceLoop_TD p3,
	      set_sticky_params_TD p4,
	      set_gb_pair_params_TD p5,
	      newEAMtype_TD p6 );

  void wrapSimInfo(setFortranSim_TD p1,
		   setFortranBox_TD p2,
		   notifyFortranCutOff_TD p3 );
  
#ifdef IS_MPI
  void wrapSimParallel( setFortranMPI_TD p1 );
#endif
}

// take care of the ForceField functions

ForceFields* currentFF;
void wrapMeFF( ForceFields* thisFF ){
  
  currentFF = thisFF;
  F90_FUNC(wrapforcefield, WRAPFORCEFIELD)( wrapFF );
}

extern "C" void wrapFF( makeAtype_TD p1,
			initFortranFF_TD p2,
			doForceLoop_TD p3,
			set_sticky_params_TD p4,
			set_gb_pair_params_TD p5,
			newEAMtype_TD p6 ){
  
  makeAtype = p1;
  initFortranFF = p2;
  currentFF->setFortranForceLoop( p3 );
  set_sticky_params = p4;
  set_gb_pair_params = p5;
  newEAMtype = p6;
}


// wrap the SimInfo functions

SimInfo* currentPlug;
void wrapMeSimInfo( SimInfo* thePlug ){
  
  currentPlug = thePlug;
  F90_FUNC(wrapsimmod, WRAPSIMMOD) ( wrapSimInfo );
}



extern "C" void wrapSimInfo( setFortranSim_TD p1,
			     setFortranBox_TD p2,
			     notifyFortranCutOff_TD p3 ){
  
  currentPlug->setInternal( p1, p2, p3 );
}



#ifdef IS_MPI

// wrap the mpiSim functions

mpiSimulation* currentMPIsim;
void wrapMeSimParallel( mpiSimulation* thisMPIsim ){
  
  currentMPIsim = thisMPIsim;
  F90_FUNC(wrapsimparallelmod, WRAPSIMPARALLELMOD) ( wrapSimParallel );
}

extern "C" void wrapSimParallel( setFortranMPI_TD p1 ){
  
  currentMPIsim->setInternal( p1 );
}


#endif // is_mpi

