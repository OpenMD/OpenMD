#include "OOPSEMinimizer.hpp"
#include "Utility.hpp"
void PRCGMinimizer::init(){

  calcG();
  
  for(int i = 0; i < direction.size(); i++){    
    direction[i] = -curG[i];
  }

}

int PRCGMinimizer::step(){
  int lsStatus;
  
  prevF = curF;
  prevG = curG;
  prevX = curX;

  //optimize along the search direction and reset minimum point value
    lsStatus = doLineSearch(direction, stepSize);

  if (lsStatus < 0)
    return -1;
  else
    return 1;
}

void PRCGMinimizer::prepareStep(){
  vector<double> deltaGrad;
  double beta;
  size_t i;

  deltaGrad.resize(ndim);
    
  //calculate the new direction using Polak-Ribiere Conjugate Gradient
  
  for(i = 0; i < curG.size(); i++)
    deltaGrad[i] = curG[i] - prevG[i];

#ifndef IS_MPI
  beta = dotProduct(deltaGrad, curG) / dotProduct(prevG, prevG);
#else
  double localDP1;
  double localDP2;
  double globalDP1;
  double globalDP2;

  localDP1 =  dotProduct(deltaGrad, curG);
  localDP2 = dotProduct(prevG, prevG);

  MPI_Allreduce(&localDP1, &globalDP1, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&localDP2, &globalDP2, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
  
  beta = globalDP1 / globalDP2;
#endif

  for(i = 0; i < direction.size(); i++)  
    direction[i] = -curG[i] + beta * direction[i];

}
