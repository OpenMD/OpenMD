#include <math.h>
#include "OOPSEMinimizer.hpp"
#include "Integrator.cpp"

OOPSEMinimizer::OOPSEMinimizer( SimInfo *theInfo, ForceFields* the_ff ,
                                              MinimizerParameterSet * param) :
              RealIntegrator(theInfo, the_ff), bShake(true), bVerbose(false) {
  dumpOut = NULL;
  statOut = NULL;

  tStats = new Thermo(info);

  
  paramSet = param;

  calcDim();
  
  curX = getCoor();
  curG.resize(ndim);

  preMove();
}

OOPSEMinimizer::~OOPSEMinimizer(){
  delete tStats;
  if(dumpOut)
    delete dumpOut;
  if(statOut)
    delete statOut;
  delete paramSet;
}

void OOPSEMinimizer::calcEnergyGradient(vector<double>& x, vector<double>& grad,
                                                                    double& energy, int& status){
 
  DirectionalAtom* dAtom;
  int index;
  double force[3];
  double dAtomGrad[6];
  int shakeStatus;

  status = 1;
  
  setCoor(x);

  if (nConstrained && bShake){
    shakeStatus = shakeR();
  }

  calcForce(1, 1);

  if (nConstrained && bShake){
    shakeStatus = shakeF();
  }

  x = getCoor();
  

  index = 0;

  for(int i = 0; i < integrableObjects.size(); i++){

    if (integrableObjects[i]->isDirectional()) {

      integrableObjects[i]->getGrad(dAtomGrad);

      //gradient is equal to -f
      grad[index++] = -dAtomGrad[0];
      grad[index++] = -dAtomGrad[1];
      grad[index++] = -dAtomGrad[2];
      grad[index++] = -dAtomGrad[3];
      grad[index++] = -dAtomGrad[4];
      grad[index++] = -dAtomGrad[5];

    }
    else{
      integrableObjects[i]->getFrc(force);

      grad[index++] = -force[0];
      grad[index++] = -force[1];
      grad[index++] = -force[2];

    }
    
  }
 
  energy = tStats->getPotential();

}

void OOPSEMinimizer::setCoor(vector<double>& x){

  DirectionalAtom* dAtom;
  int index;
  double position[3];
  double eulerAngle[3];


  index = 0;
  
  for(int i = 0; i < integrableObjects.size(); i++){
    
    position[0] = x[index++];
    position[1] = x[index++];
    position[2] = x[index++];

    integrableObjects[i]->setPos(position);

    if (integrableObjects[i]->isDirectional()){
 
      eulerAngle[0] = x[index++];
      eulerAngle[1] = x[index++];
      eulerAngle[2] = x[index++];

      integrableObjects[i]->setEuler(eulerAngle[0], 
                                     eulerAngle[1], 
                                     eulerAngle[2]);

    }
    
  }
  
}

vector<double> OOPSEMinimizer::getCoor(){
 
  DirectionalAtom* dAtom;
  int index;
  double position[3];
  double eulerAngle[3];
  vector<double> x;

  x.resize(getDim());

  index = 0;
  
  for(int i = 0; i < integrableObjects.size(); i++){
    integrableObjects[i]->getPos(position);

    x[index++] = position[0];
    x[index++] = position[1];
    x[index++] = position[2];

    if (integrableObjects[i]->isDirectional()){

      integrableObjects[i]->getEulerAngles(eulerAngle);
      
      x[index++] = eulerAngle[0];
      x[index++] = eulerAngle[1];
      x[index++] = eulerAngle[2];
       
    }
    
  }

  return x;

}


int OOPSEMinimizer::shakeR(){
  int i, j;
  int done;
  double posA[3], posB[3];
  double velA[3], velB[3];
  double pab[3];
  double rab[3];
  int a, b, ax, ay, az, bx, by, bz;
  double rma, rmb;
  double dx, dy, dz;
  double rpab;
  double rabsq, pabsq, rpabsq;
  double diffsq;
  double gab;
  int iteration;

  for (i = 0; i < nAtoms; i++){
    moving[i] = 0;
    moved[i] = 1;
  }

  iteration = 0;
  done = 0;
  while (!done && (iteration < maxIteration)){
    done = 1;
    for (i = 0; i < nConstrained; i++){
      a = constrainedA[i];
      b = constrainedB[i];

      ax = (a * 3) + 0;
      ay = (a * 3) + 1;
      az = (a * 3) + 2;

      bx = (b * 3) + 0;
      by = (b * 3) + 1;
      bz = (b * 3) + 2;

      if (moved[a] || moved[b]){
        atoms[a]->getPos(posA);
        atoms[b]->getPos(posB);

        for (j = 0; j < 3; j++)
          pab[j] = posA[j] - posB[j];

        //periodic boundary condition

        info->wrapVector(pab);

        pabsq = pab[0] * pab[0] + pab[1] * pab[1] + pab[2] * pab[2];

        rabsq = constrainedDsqr[i];
        diffsq = rabsq - pabsq;

        // the original rattle code from alan tidesley
        if (fabs(diffsq) > (tol * rabsq * 2)){
          rab[0] = oldPos[ax] - oldPos[bx];
          rab[1] = oldPos[ay] - oldPos[by];
          rab[2] = oldPos[az] - oldPos[bz];

          info->wrapVector(rab);
          
          rpab = rab[0] * pab[0] + rab[1] * pab[1] + rab[2] * pab[2];

          rpabsq = rpab * rpab;


          if (rpabsq < (rabsq * -diffsq)){
#ifdef IS_MPI
            a = atoms[a]->getGlobalIndex();
            b = atoms[b]->getGlobalIndex();
#endif //is_mpi
            //cerr << "Waring: constraint failure" << endl;
            gab = sqrt(rabsq/pabsq);
            rab[0] = (posA[0] - posB[0])*gab;
            rab[1]= (posA[1] - posB[1])*gab;
            rab[2] = (posA[2] - posB[2])*gab;
            
            info->wrapVector(rab);
            
            rpab = rab[0] * pab[0] + rab[1] * pab[1] + rab[2] * pab[2];
            
          }

          //rma = 1.0 / atoms[a]->getMass();
          //rmb = 1.0 / atoms[b]->getMass();
          rma = 1.0;
          rmb =1.0;

          gab = diffsq / (2.0 * (rma + rmb) * rpab);

          dx = rab[0] * gab;
          dy = rab[1] * gab;
          dz = rab[2] * gab;

          posA[0] += rma * dx;
          posA[1] += rma * dy;
          posA[2] += rma * dz;

          atoms[a]->setPos(posA);

          posB[0] -= rmb * dx;
          posB[1] -= rmb * dy;
          posB[2] -= rmb * dz;

          atoms[b]->setPos(posB);

          moving[a] = 1;
          moving[b] = 1;
          done = 0;
        }
      }
    }

    for (i = 0; i < nAtoms; i++){
      moved[i] = moving[i];
      moving[i] = 0;
    }

    iteration++;
  }

  if (!done){
    cerr << "Waring: can not constraint within maxIteration" << endl;
    return -1;
  }
  else
    return 1;
}


//remove constraint force along the bond direction
int OOPSEMinimizer::shakeF(){
  int i, j;
  int done;
  double posA[3], posB[3];
  double frcA[3], frcB[3];
  double rab[3], fpab[3];
  int a, b, ax, ay, az, bx, by, bz;
  double rma, rmb;
  double rvab;
  double gab;
  double rabsq;
  double rfab;
  int iteration;

  for (i = 0; i < nAtoms; i++){
    moving[i] = 0;
    moved[i] = 1;
  }

  done = 0;
  iteration = 0;
  while (!done && (iteration < maxIteration)){
    done = 1;

    for (i = 0; i < nConstrained; i++){
      a = constrainedA[i];
      b = constrainedB[i];

      ax = (a * 3) + 0;
      ay = (a * 3) + 1;
      az = (a * 3) + 2;

      bx = (b * 3) + 0;
      by = (b * 3) + 1;
      bz = (b * 3) + 2;

      if (moved[a] || moved[b]){

        atoms[a]->getPos(posA);
        atoms[b]->getPos(posB);

        for (j = 0; j < 3; j++)
          rab[j] = posA[j] - posB[j];

        info->wrapVector(rab);

        atoms[a]->getFrc(frcA);
        atoms[b]->getFrc(frcB);

        //rma = 1.0 / atoms[a]->getMass();
        //rmb = 1.0 / atoms[b]->getMass();
        rma = 1.0;
        rmb = 1.0;
        
        
        fpab[0] = frcA[0] * rma - frcB[0] * rmb;
        fpab[1] = frcA[1] * rma - frcB[1] * rmb;
        fpab[2] = frcA[2] * rma - frcB[2] * rmb;


          gab=fpab[0] * fpab[0] + fpab[1] * fpab[1] + fpab[2] * fpab[2];
          
          if (gab < 1.0)
            gab = 1.0;
         
          rabsq = rab[0] * rab[0] + rab[1] * rab[1] + rab[2] * rab[2];
          rfab = rab[0] * fpab[0] + rab[1] * fpab[1] + rab[2] * fpab[2];

          if (fabs(rfab) > sqrt(rabsq*gab) * 0.00001){

            gab = -rfab /(rabsq*(rma + rmb));
            
            frcA[0] = rab[0] * gab;
            frcA[1] = rab[1] * gab;
            frcA[2] = rab[2] * gab;

            atoms[a]->addFrc(frcA);
            

            frcB[0] = -rab[0] * gab;
            frcB[1] = -rab[1] * gab;
            frcB[2] = -rab[2] * gab;

            atoms[b]->addFrc(frcB);
          
            moving[a] = 1;
            moving[b] = 1;
            done = 0;
          }             
      }
    }

    for (i = 0; i < nAtoms; i++){
      moved[i] = moving[i];
      moving[i] = 0;
    }

    iteration++;
  }

  if (!done){
    cerr << "Waring: can not constraint within maxIteration" << endl;
    return -1;
  }
  else
    return 1;
}



//calculate the value of object function
void OOPSEMinimizer::calcF(){
  calcEnergyGradient(curX, curG, curF, egEvalStatus);
}
    
void OOPSEMinimizer::calcF(vector<double>& x, double&f, int& status){
  vector<double> tempG;
  tempG.resize(x.size());
  
  calcEnergyGradient(x, tempG, f, status);
}

//calculate the gradient
void OOPSEMinimizer::calcG(){
  calcEnergyGradient(curX, curG, curF, egEvalStatus);
}

void OOPSEMinimizer::calcG(vector<double>& x, vector<double>& g, double& f, int& status){
  calcEnergyGradient(x, g, f, status);
}

void OOPSEMinimizer::calcDim(){
  DirectionalAtom* dAtom;

  ndim = 0;

  for(int i = 0; i < integrableObjects.size(); i++){
    ndim += 3;
    if (integrableObjects[i]->isDirectional())
      ndim += 3;      
  }
}

void OOPSEMinimizer::setX(vector < double > & x){

    if (x.size() != ndim && bVerbose){
      //sprintf(painCave.errMsg,
      //          "OOPSEMinimizer Error: dimesion of x and curX does not match\n");
     // painCave.isFatal = 1;
     // simError();
    }

    curX = x;
}

void OOPSEMinimizer::setG(vector < double > & g){

    if (g.size() != ndim && bVerbose){
      //sprintf(painCave.errMsg,
      //          "OOPSEMinimizer Error: dimesion of g and curG does not match\n");
     // painCave.isFatal = 1;
      //simError();
    }

    curG = g;
}

void OOPSEMinimizer::writeOut(vector<double>& x, double iter){

  setX(x);

  calcG();

  dumpOut->writeDump(iter);
  statOut->writeStat(iter);
}


void OOPSEMinimizer::printMinimizerInfo(){
  cout << "--------------------------------------------------------------------" << endl;
  cout << minimizerName << endl;
  cout << "minimization parameter set" << endl;
  cout << "function tolerance = " << paramSet->getFTol() << endl;
  cout << "gradient tolerance = " << paramSet->getGTol() << endl;
  cout << "step tolerance = "<< paramSet->getFTol() << endl;
  cout << "absolute gradient tolerance = " << endl;
  cout << "max iteration = " << paramSet->getMaxIteration() << endl;
  cout << "max line search iteration = " << paramSet->getLineSearchMaxIteration() <<endl;
  cout << "shake algorithm = " << bShake << endl;
  cout << "--------------------------------------------------------------------" << endl;

}

/**
 * In thoery, we need to find the minimum along the search direction
 * However, function evaluation is too expensive. 
 * At the very begining of the problem, we check the search direction and make sure
 * it is a descent direction
 * we will compare the energy of two end points,
 * if the right end point has lower energy, we just take it
 *
 *
 *
 */

int OOPSEMinimizer::doLineSearch(vector<double>& direction, double stepSize){
  vector<double> xa;
  vector<double> xb;
  vector<double> xc;
  vector<double> ga;
  vector<double> gb;
  vector<double> gc;
  double fa;
  double fb;
  double fc;
  double a;
  double b;
  double c;
  int status;
  double initSlope;
  double slopeA;
  double slopeB;
  double slopeC;
  bool foundLower;
  int iter;
  int maxLSIter;
  double mu;
  double eta;
  double ftol;  
  double lsTol;

  xa.resize(ndim);
  xb.resize(ndim);
  xc.resize(ndim);

  ga.resize(ndim);
  gb.resize(ndim);
  gc.resize(ndim);

  a = 0.0;
  fa =  curF;    
  xa = curX;
  ga = curG;
  c = a + stepSize; 
  ftol = paramSet->getFTol();
  lsTol = paramSet->getLineSearchTol();
        
  //calculate the derivative at a = 0 
  slopeA = 0;
  for (size_t i = 0; i < ndim; i++)
    slopeA += curG[i]*direction[i];

  initSlope = slopeA;
  
  // if  going uphill, use negative gradient as searching direction 
  if (slopeA > 0) {

    if (bVerbose){
      cout << "LineSearch Warning: initial searching direction is not a descent searching direction, "
             << " use negative gradient instead. Therefore, finding a smaller vaule of function "
             << " is guaranteed"
             << endl;
    }    
    
    for (size_t i = 0; i < ndim; i++)
      direction[i] = -curG[i];    
    
    for (size_t i = 0; i < ndim; i++)
      slopeA += curG[i]*direction[i];
    
    initSlope = slopeA;
  }
    
  // Take a trial step
  for(size_t i = 0; i < ndim; i++)
    xc[i] = curX[i] + direction[i] * c;      

  calcG(xc, gc, fc, status);

  if (status < 0){
    if (bVerbose)
      cerr << "Function Evaluation Error" << endl;
  }

  //calculate the derivative at c
  slopeC = 0;
  for (size_t i = 0; i < ndim; i++)
    slopeC += gc[i]*direction[i];

  // found a lower point
  if (fc < fa) {
    curX = xc;
    curG = gc;
    curF = fc;
    return LS_SUCCEED;
  }
  else {

    if (slopeC > 0)
    stepSize *= 0.618034;
  }    

  maxLSIter = paramSet->getLineSearchMaxIteration();
   
  iter = 0;
  
  do {
    // Select a new trial point.
    // If the derivatives at points a & c have different sign we use cubic interpolate    
    //if (slopeC > 0){     
      eta = 3 *(fa -fc) /(c - a) + slopeA + slopeC;
      mu = sqrt(eta * eta - slopeA * slopeC);      
      b = a + (c - a) * (1 - (slopeC + mu - eta) /(slopeC - slopeA + 2 * mu));       

      if (b < lsTol){
        if (bVerbose)
          cout << "stepSize is less than line search tolerance" << endl;
        break;        
      }
    //}

    // Take a trial step to this new point - new coords in xb 
    for(size_t i = 0; i < ndim; i++)
      xb[i] = curX[i] + direction[i] * b;      

    //function evaluation
    calcG(xb, gb, fb, status);

    if (status < 0){
      if (bVerbose)
        cerr << "Function Evaluation Error" << endl;
    }

  //calculate the derivative at c
    slopeB = 0;
    for (size_t i = 0; i < ndim; i++)
      slopeB += gb[i]*direction[i];

    //Amijo Rule to stop the line search 
    if (fb <= curF +  initSlope * ftol * b) {
      curF = fb;
      curX = xb;
      curG = gb;
      return LS_SUCCEED;
     }
 
    if (slopeB <0 &&  fb < fa) {
      //replace a by b
      fa = fb;
      a = b;
      slopeA = slopeB;

      // swap coord  a/b 
      std::swap(xa, xb);
      std::swap(ga, gb);
    }
    else {
      //replace c by b
      fc = fb;
      c = b;
      slopeC = slopeB;

      // swap coord  b/c 
      std::swap(gb, gc);
      std::swap(xb, xc);
    }
     

     iter++;
  } while((fb > fa || fb > fc) && (iter < maxLSIter));     
  
  if(fb < curF || iter >= maxLSIter) {
    //could not find a lower value, we might just go uphill.      
    return LS_ERROR;
  }

  //select the end point

  if (fa <= fc) {
    curX = xa;
    curG = ga;
    curF = fa;
  }
  else {
    curX = xc;
    curG = gc;
    curF = fc;    
  }

  return LS_SUCCEED;
    
}

void OOPSEMinimizer::minimize(){

  int convgStatus;
  int stepStatus;
  int maxIter;
  int writeFrq;
  int nextWriteIter;
  
  if (bVerbose)
    printMinimizerInfo();

  dumpOut = new DumpWriter(info);
  statOut = new StatWriter(info);
  
  init();

  writeFrq = paramSet->getWriteFrq();
  nextWriteIter = writeFrq;
  
  maxIter = paramSet->getMaxIteration();
  
  for (curIter = 1; curIter <= maxIter; curIter++){

    stepStatus = step();

    if (nConstrained && bShake)
      preMove();
    
    if (stepStatus < 0){
      saveResult();
      minStatus = MIN_LSERROR;
      cerr << "OOPSEMinimizer Error: line search error, please try a small stepsize" << endl;
      return;
    }

    if (curIter == nextWriteIter){
      nextWriteIter += writeFrq;
      writeOut(curX, curIter); 
    }

    convgStatus = checkConvg();

    if (convgStatus > 0){
      saveResult();
      minStatus = MIN_CONVERGE;
      return;
    }
    
    prepareStep();

  }


  if (bVerbose) {
    cout << "OOPSEMinimizer Warning: "
           << minimizerName << " algorithm did not converge within "
           << maxIter << " iteration" << endl;
  }
  minStatus = MIN_MAXITER;
  saveResult();
  
}
