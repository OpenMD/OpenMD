#include <math.h>
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

// Basic isotropic thermostating and barostating via the Melchionna
// modification of the Hoover algorithm:
//
//    Melchionna, S., Ciccotti, G., and Holian, B. L., 1993,
//       Molec. Phys., 78, 533.
//
//           and
//
//    Hoover, W. G., 1986, Phys. Rev. A, 34, 2499.

template<typename T> NPTi<T>::NPTi ( SimInfo *theInfo, ForceFields* the_ff):
  T( theInfo, the_ff )
{
  GenericData* data;
  DoubleVectorGenericData * etaValue;
  vector<double> etaArray;

  eta = 0.0;
  oldEta = 0.0;

  if( theInfo->useInitXSstate ){
    // retrieve eta from simInfo if
    data = info->getProperty(ETAVALUE_ID);
    if(data){
      etaValue = dynamic_cast<DoubleVectorGenericData*>(data);
      
      if(etaValue){
	eta = (*etaValue)[0];
	oldEta = eta;
      }
    }
  }
}

template<typename T> NPTi<T>::~NPTi() {
  //nothing for now
}

template<typename T> void NPTi<T>::resetIntegrator() {
  eta = 0.0;
  T::resetIntegrator();
}

template<typename T> void NPTi<T>::evolveEtaA() {
  eta += dt2 * ( instaVol * (instaPress - targetPressure) /
		 (p_convert*NkBT*tb2));
  oldEta = eta;
}

template<typename T> void NPTi<T>::evolveEtaB() {

  prevEta = eta;
  eta = oldEta + dt2 * ( instaVol * (instaPress - targetPressure) /
		 (p_convert*NkBT*tb2));
}

template<typename T> void NPTi<T>::calcVelScale(void) {
  vScale = chi + eta;
}

template<typename T> void NPTi<T>::getVelScaleA(double sc[3], double vel[3]) {
  int i;

  for(i=0; i<3; i++) sc[i] = vel[i] * vScale;
}

template<typename T> void NPTi<T>::getVelScaleB(double sc[3], int index ){
  int i;

  for(i=0; i<3; i++) sc[i] = oldVel[index*3 + i] * vScale;
}


template<typename T> void NPTi<T>::getPosScale(double pos[3], double COM[3],
					       int index, double sc[3]){
  int j;

  for(j=0; j<3; j++)
    sc[j] = ( oldPos[index*3+j] + pos[j]) / 2.0 - COM[j];

  for(j=0; j<3; j++)
    sc[j] *= eta;
}

template<typename T> void NPTi<T>::scaleSimBox( void ){

  double scaleFactor;

  scaleFactor = exp(dt*eta);

  if ((scaleFactor > 1.1) || (scaleFactor < 0.9)) {
    sprintf( painCave.errMsg,
             "NPTi error: Attempting a Box scaling of more than 10 percent"
             " check your tauBarostat, as it is probably too small!\n"
             " eta = %lf, scaleFactor = %lf\n", eta, scaleFactor
             );
    painCave.isFatal = 1;
    simError();
  } else {
    info->scaleBox(scaleFactor);
  }

}

template<typename T> bool NPTi<T>::etaConverged() {

  return ( fabs(prevEta - eta) <= etaTolerance );
}

template<typename T> double NPTi<T>::getConservedQuantity(void){

  double conservedQuantity;
  double Energy;
  double thermostat_kinetic;
  double thermostat_potential;
  double barostat_kinetic;
  double barostat_potential;

  Energy = tStats->getTotalE();

  thermostat_kinetic = fkBT* tt2 * chi * chi /
    (2.0 * eConvert);

  thermostat_potential = fkBT* integralOfChidt / eConvert;


  barostat_kinetic = 3.0 * NkBT * tb2 * eta * eta /
    (2.0 * eConvert);

  barostat_potential = (targetPressure * tStats->getVolume() / p_convert) /
    eConvert;

  conservedQuantity = Energy + thermostat_kinetic + thermostat_potential +
    barostat_kinetic + barostat_potential;

//   cout.width(8);
//   cout.precision(8);

//   cerr << info->getTime() << "\t" << Energy << "\t" << thermostat_kinetic <<
//       "\t" << thermostat_potential << "\t" << barostat_kinetic <<
//       "\t" << barostat_potential << "\t" << conservedQuantity << endl;
  return conservedQuantity;
}

template<typename T> string NPTi<T>::getAdditionalParameters(void){
  string parameters;
  const int BUFFERSIZE = 2000; // size of the read buffer
  char buffer[BUFFERSIZE];

  sprintf(buffer,"\t%G\t%G;", chi, integralOfChidt);
  parameters += buffer;

  sprintf(buffer,"\t%G\t0\t0;", eta);
  parameters += buffer;

  sprintf(buffer,"\t0\t%G\t0;", eta);
  parameters += buffer;

  sprintf(buffer,"\t0\t0\t%G;", eta);
  parameters += buffer;

  return parameters;

}
