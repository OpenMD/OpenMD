#include "OOPSEMinimizer.hpp"
#include "Utility.hpp"

SDMinimizer::SDMinimizer(SimInfo *theInfo, ForceFields* the_ff ,
                                                               MinimizerParameterSet * param)
                         :OOPSEMinimizer(theInfo, the_ff, param){

  direction.resize(ndim);
  stepSize = paramSet->getStepSize();                         
}

void SDMinimizer::init(){

  calcG();
  
  for(int i = 0; i < direction.size(); i++){    
    direction[i] = -curG[i];
  }  
}  

int SDMinimizer::step(){
  int lsStatus;
 
  prevF = curF;
  
  //optimize along the search direction and reset minimum point value
    lsStatus = doLineSearch(direction, stepSize);

  if (lsStatus < 0)
    return -1;
  else
    return 1;
}

void SDMinimizer::prepareStep(){
  
  for(int i = 0; i < direction.size(); i++){    
    direction[i] = -curG[i];
  }  
}  
int SDMinimizer::checkConvg(){
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

  return CONVG_UNCONVG;
}

