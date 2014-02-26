#ifdef IS_MPI
#include <mpi.h>
#endif

#include "config.h"
#include "optimization/Problem.hpp"

namespace QuantLib {
  RealType Problem::DotProduct(DynamicVector<RealType>& v1, 
                               DynamicVector<RealType>& v2){ 
    RealType dp = dot(v1, v2);
#ifdef IS_MPI
    // in parallel, we need to add up the contributions from all
    // processors:
    MPI_Allreduce(MPI_IN_PLACE, &dp, 1, MPI_REALTYPE, 
                  MPI_SUM, MPI_COMM_WORLD);
#endif
    return dp;    
  }

  RealType Problem::computeGradientNormValue(DynamicVector<RealType>& grad_f) { 
    
    RealType dot = grad_f.lengthSquare();
#ifdef IS_MPI
    // in parallel, we need to add up the contributions from all
    // processors:
    MPI_Allreduce(MPI_IN_PLACE, &dot, 1, MPI_REALTYPE, 
                  MPI_SUM, MPI_COMM_WORLD);
#endif
    return dot;

  } 
}
