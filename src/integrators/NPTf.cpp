#include <math.h>

#include "math/MatVec3.h"
#include "primitives/Atom.hpp"
#include "primitives/SRI.hpp"
#include "primitives/AbstractClasses.hpp"
#include "brains/SimInfo.hpp"
#include "UseTheForce/ForceFields.hpp"
#include "brains/Thermo.hpp"
#include "io/ReadWrite.hpp"
#include "integrators/Integrator.hpp"
#include "utils/simError.h"

#ifdef IS_MPI
#include "brains/mpiSimulation.hpp"
#endif

// Basic non-isotropic thermostating and barostating via the Melchionna
// modification of the Hoover algorithm:
//
//    Melchionna, S., Ciccotti, G., and Holian, B. L., 1993,
//       Molec. Phys., 78, 533.
//
//           and
//
//    Hoover, W. G., 1986, Phys. Rev. A, 34, 2499.

template<typename T> NPTf<T>::NPTf ( SimInfo *theInfo, ForceFields* the_ff):
  T( theInfo, the_ff )
{
  GenericData* data;
  DoubleArrayData * etaValue;
  vector<double> etaArray;
  int i,j;

  for(i = 0; i < 3; i++){
    for (j = 0; j < 3; j++){

      eta[i][j] = 0.0;
      oldEta[i][j] = 0.0;
    }
  }


  if( theInfo->useInitXSstate ){
    // retrieve eta array from simInfo if it exists
    data = info->getProperty(ETAVALUE_ID);
    if(data){
      etaValue = dynamic_cast<DoubleArrayData*>(data);
      
      if(etaValue){
	etaArray = etaValue->getData();
	
	for(i = 0; i < 3; i++){
	  for (j = 0; j < 3; j++){
	    eta[i][j] = etaArray[3*i+j];
	    oldEta[i][j] = eta[i][j];
	  }
	}
      }
    }
  }

}

template<typename T> NPTf<T>::~NPTf() {

  // empty for now
}

template<typename T> void NPTf<T>::resetIntegrator() {

  int i, j;

  for(i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      eta[i][j] = 0.0;

  T::resetIntegrator();
}

template<typename T> void NPTf<T>::evolveEtaA() {

  int i, j;

  for(i = 0; i < 3; i ++){
    for(j = 0; j < 3; j++){
      if( i == j)
        eta[i][j] += dt2 *  instaVol *
	  (press[i][j] - targetPressure/p_convert) / (NkBT*tb2);
      else
        eta[i][j] += dt2 * instaVol * press[i][j] / (NkBT*tb2);
    }
  }
  
  for(i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      oldEta[i][j] = eta[i][j];
}

template<typename T> void NPTf<T>::evolveEtaB() {

  int i,j;

  for(i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      prevEta[i][j] = eta[i][j];

  for(i = 0; i < 3; i ++){
    for(j = 0; j < 3; j++){
      if( i == j) {
	eta[i][j] = oldEta[i][j] + dt2 *  instaVol *
	  (press[i][j] - targetPressure/p_convert) / (NkBT*tb2);
      } else {
	eta[i][j] = oldEta[i][j] + dt2 * instaVol * press[i][j] / (NkBT*tb2);
      }
    }
  }
}

template<typename T> void NPTf<T>::calcVelScale(void){
  int i,j;

  for (i = 0; i < 3; i++ ) {
    for (j = 0; j < 3; j++ ) {
      vScale[i][j] = eta[i][j];

      if (i == j) {
        vScale[i][j] += chi;
      }
    }
  }
}

template<typename T> void NPTf<T>::getVelScaleA(double sc[3], double vel[3]) {
 
  matVecMul3( vScale, vel, sc );
}

template<typename T> void NPTf<T>::getVelScaleB(double sc[3], int index ){
  int j;
  double myVel[3];

  for (j = 0; j < 3; j++)
    myVel[j] = oldVel[3*index + j];
  
  matVecMul3( vScale, myVel, sc );
}

template<typename T> void NPTf<T>::getPosScale(double pos[3], double COM[3],
					       int index, double sc[3]){
  int j;
  double rj[3];

  for(j=0; j<3; j++)
    rj[j] = ( oldPos[index*3+j] + pos[j]) / 2.0 - COM[j];

  matVecMul3( eta, rj, sc );
}

template<typename T> void NPTf<T>::scaleSimBox( void ){

  int i,j,k;
  double scaleMat[3][3];
  double eta2ij;
  double bigScale, smallScale, offDiagMax;
  double hm[3][3], hmnew[3][3];



  // Scale the box after all the positions have been moved:

  // Use a taylor expansion for eta products:  Hmat = Hmat . exp(dt * etaMat)
  //  Hmat = Hmat . ( Ident + dt * etaMat  + dt^2 * etaMat*etaMat / 2)

  bigScale = 1.0;
  smallScale = 1.0;
  offDiagMax = 0.0;

  for(i=0; i<3; i++){
    for(j=0; j<3; j++){

      // Calculate the matrix Product of the eta array (we only need
      // the ij element right now):

      eta2ij = 0.0;
      for(k=0; k<3; k++){
        eta2ij += eta[i][k] * eta[k][j];
      }

      scaleMat[i][j] = 0.0;
      // identity matrix (see above):
      if (i == j) scaleMat[i][j] = 1.0;
      // Taylor expansion for the exponential truncated at second order:
      scaleMat[i][j] += dt*eta[i][j]  + 0.5*dt*dt*eta2ij;
      

      if (i != j)
        if (fabs(scaleMat[i][j]) > offDiagMax)
          offDiagMax = fabs(scaleMat[i][j]);
    }

    if (scaleMat[i][i] > bigScale) bigScale = scaleMat[i][i];
    if (scaleMat[i][i] < smallScale) smallScale = scaleMat[i][i];
  }

  if ((bigScale > 1.01) || (smallScale < 0.99)) {
    sprintf( painCave.errMsg,
             "NPTf error: Attempting a Box scaling of more than 1 percent.\n"
             " Check your tauBarostat, as it is probably too small!\n\n"
             " scaleMat = [%lf\t%lf\t%lf]\n"
             "            [%lf\t%lf\t%lf]\n"
             "            [%lf\t%lf\t%lf]\n"
             "      eta = [%lf\t%lf\t%lf]\n"
             "            [%lf\t%lf\t%lf]\n"
             "            [%lf\t%lf\t%lf]\n",
             scaleMat[0][0],scaleMat[0][1],scaleMat[0][2],
             scaleMat[1][0],scaleMat[1][1],scaleMat[1][2],
             scaleMat[2][0],scaleMat[2][1],scaleMat[2][2],
             eta[0][0],eta[0][1],eta[0][2],
             eta[1][0],eta[1][1],eta[1][2],
             eta[2][0],eta[2][1],eta[2][2]);
    painCave.isFatal = 1;
    simError();
  } else if (offDiagMax > 0.01) {
    sprintf( painCave.errMsg,
             "NPTf error: Attempting an off-diagonal Box scaling of more than 1 percent.\n"
             " Check your tauBarostat, as it is probably too small!\n\n"
             " scaleMat = [%lf\t%lf\t%lf]\n"
             "            [%lf\t%lf\t%lf]\n"
             "            [%lf\t%lf\t%lf]\n"
             "      eta = [%lf\t%lf\t%lf]\n"
             "            [%lf\t%lf\t%lf]\n"
             "            [%lf\t%lf\t%lf]\n",
             scaleMat[0][0],scaleMat[0][1],scaleMat[0][2],
             scaleMat[1][0],scaleMat[1][1],scaleMat[1][2],
             scaleMat[2][0],scaleMat[2][1],scaleMat[2][2],
             eta[0][0],eta[0][1],eta[0][2],
             eta[1][0],eta[1][1],eta[1][2],
             eta[2][0],eta[2][1],eta[2][2]);
    painCave.isFatal = 1;
    simError();
  } else {
    info->getBoxM(hm);
    matMul3(hm, scaleMat, hmnew);
    info->setBoxM(hmnew);
  }
}

template<typename T> bool NPTf<T>::etaConverged() {
  int i;
  double diffEta, sumEta;

  sumEta = 0;
  for(i = 0; i < 3; i++)
    sumEta += pow(prevEta[i][i] - eta[i][i], 2);

  diffEta = sqrt( sumEta / 3.0 );

  return ( diffEta <= etaTolerance );
}

template<typename T> double NPTf<T>::getConservedQuantity(void){

  double conservedQuantity;
  double totalEnergy;
  double thermostat_kinetic;
  double thermostat_potential;
  double barostat_kinetic;
  double barostat_potential;
  double trEta;
  double a[3][3], b[3][3];

  totalEnergy = tStats->getTotalE();

  thermostat_kinetic = fkBT * tt2 * chi * chi /
    (2.0 * eConvert);

  thermostat_potential = fkBT* integralOfChidt / eConvert;

  transposeMat3(eta, a);
  matMul3(a, eta, b);
  trEta = matTrace3(b);

  barostat_kinetic = NkBT * tb2 * trEta /
    (2.0 * eConvert);

  barostat_potential = (targetPressure * tStats->getVolume() / p_convert) /
    eConvert;

  conservedQuantity = totalEnergy + thermostat_kinetic + thermostat_potential +
    barostat_kinetic + barostat_potential;

  return conservedQuantity;

}

template<typename T> string NPTf<T>::getAdditionalParameters(void){
  string parameters;
  const int BUFFERSIZE = 2000; // size of the read buffer
  char buffer[BUFFERSIZE];

  sprintf(buffer,"\t%G\t%G;", chi, integralOfChidt);
  parameters += buffer;

  for(int i = 0; i < 3; i++){
    sprintf(buffer,"\t%G\t%G\t%G;", eta[i][0], eta[i][1], eta[i][2]);
    parameters += buffer;
  }

  return parameters;

}
