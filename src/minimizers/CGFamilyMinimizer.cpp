#include <math.h>

#include "OOPSEMinimizer.hpp"
#include "Utility.hpp"

CGFamilyMinmizer::CGFamilyMinmizer(SimInfo *theInfo, ForceFields* the_ff ,
                                                               MinimizerParameterSet * param)
                         :OOPSEMinimizer(theInfo, the_ff, param){
  prevG.resize(ndim);
  prevX.resize(ndim);
  direction.resize(ndim);

  stepSize = paramSet->getStepSize();

}
int CGFamilyMinmizer::checkConvg(){
  double fTol;
  double relativeFTol;  // relative tolerance
  double deltaF;
  double gTol;
  double relativeGTol;
  double gnorm;
  

  // test function tolerance test
  fTol =paramSet->getFTol();
  relativeFTol = fTol * std::max(1.0,fabs(curF));  // relative tolerance
  deltaF = prevF - curF;
  
  if (fabs(deltaF) <= relativeFTol) {

    if (bVerbose){
      cout << "function value tolerance test passed" << endl;
      cout << "ftol = " << fTol
             << "\tdeltaf = " << deltaF<< endl;
    }
    return CONVG_FTOL;
  }
  
//gradient tolerance test
  gTol = paramSet->getGTol();
  relativeGTol = gTol * std::max(1.0,fabs(curF));

#ifndef IS_MPI
  gnorm = sqrt(dotProduct(curG, curG));
#else
  double localDP;
  double globalDP;

  localDP = dotProduct(curG, curG);
  MPI_Allreduce(&localDP, &globalDP, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);  
  gnorm  = sqrt(globalDP);
#endif

  if (gnorm <= relativeGTol) {
      cout << "gradient tolerance test" << endl;
      cout << "gnorm = " << gnorm
             << "\trelativeGTol = " << relativeGTol<< endl;
    return CONVG_GTOL;
  }
  
  //absolute gradient tolerance test

  if (gnorm <= gTol) {
      cout << "absolute gradient tolerance test" << endl;
      cout << "gnorm = " << gnorm
             << "\tgTol = " << gTol<< endl;
    return CONVG_ABSGTOL;
  }

/*
  //test step tolerance test

  double stepTol = paramSet->getStepTol();
  double snorm = stepTolNorm();
  double xnorm =  Norm2(curX);
  double sTol  = stepTol*max(1.0, xnorm);
  if (snorm  <= sTol) {
      cout << "step tolerance test passed" << endl;
      cout << "stol = " << sTol
             << "\tsnorm = " << snorm<< endl;
    return CONVG_STEPTOL;
  }
*/   
  // did not converge yet
  return CONVG_UNCONVG;
}

