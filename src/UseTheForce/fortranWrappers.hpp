#ifndef __FORTRANWRAPPERS_H__
#define __FORTRANWRAPPERS_H__

#include "fortranWrapDefines.hpp"
#include "ForceFields.hpp"
#include "SimInfo.hpp"

#ifdef IS_MPI
#include "mpiSimulation.hpp"
#endif // is_mpi


extern makeAtype_TD          makeAtype;
extern newEAMtype_TD         newEAMtype;
extern initFortranFF_TD      initFortranFF;
extern set_sticky_params_TD  set_sticky_params;
extern set_gb_pair_params_TD set_gb_pair_params;

extern void wrapMeFF( ForceFields* thisFF );

extern void wrapMeSimInfo( SimInfo* thePlug );

#ifdef IS_MPI
extern void wrapMeSimParallel( mpiSimulation* thisMPIsim );
#endif // is_mpi


#endif // fortranWrappers.hpp
