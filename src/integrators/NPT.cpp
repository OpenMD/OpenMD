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

template<typename T> NPT<T>::NPT ( SimInfo *theInfo, ForceFields* the_ff):
  T( theInfo, the_ff )
{
  GenericData* data;
  DoubleGenericData * chiValue;
  DoubleGenericData * integralOfChidtValue;

  chiValue = NULL;
  integralOfChidtValue = NULL;

  chi = 0.0;
  integralOfChidt = 0.0;
  have_tau_thermostat = 0;
  have_tau_barostat = 0;
  have_target_temp = 0;
  have_target_pressure = 0;
  have_chi_tolerance = 0;
  have_eta_tolerance = 0;
  have_pos_iter_tolerance = 0;

  // retrieve chi and integralOfChidt from simInfo
  data = info->getProperty(CHIVALUE_ID);
  if(data){
    chiValue = dynamic_cast<DoubleGenericData*>(data);
  }

  data = info->getProperty(INTEGRALOFCHIDT_ID);
  if(data){
    integralOfChidtValue = dynamic_cast<DoubleGenericData*>(data);
  }

  // chi and integralOfChidt should appear by pair
  if(chiValue && integralOfChidtValue){
    chi = chiValue->getData();
    integralOfChidt = integralOfChidtValue->getData();
  }

  oldPos = new double[3*integrableObjects.size()];
  oldVel = new double[3*integrableObjects.size()];
  oldJi = new double[3*integrableObjects.size()];

}

template<typename T> NPT<T>::~NPT() {
  delete[] oldPos;
  delete[] oldVel;
  delete[] oldJi;
}

template<typename T> void NPT<T>::moveA() {

  //new version of NPT
  int i, j, k;
  double Tb[3], ji[3];
  double mass;
  double vel[3], pos[3], frc[3];
  double sc[3];
  double COM[3];

  instaTemp = tStats->getTemperature();
  tStats->getPressureTensor( press );
  instaPress = p_convert * (press[0][0] + press[1][1] + press[2][2]) / 3.0;
  instaVol = tStats->getVolume();

  tStats->getCOM(COM);

  //evolve velocity half step

  calcVelScale();
  
  for( i=0; i<integrableObjects.size(); i++ ){

    integrableObjects[i]->getVel( vel );
    integrableObjects[i]->getFrc( frc );

    mass = integrableObjects[i]->getMass();

    getVelScaleA( sc, vel );

    for (j=0; j < 3; j++) {

      // velocity half step  (use chi from previous step here):
      vel[j] += dt2 * ((frc[j] / mass ) * eConvert - sc[j]);

    }

    integrableObjects[i]->setVel( vel );

    if( integrableObjects[i]->isDirectional() ){

      // get and convert the torque to body frame

      integrableObjects[i]->getTrq( Tb );
      integrableObjects[i]->lab2Body( Tb );

      // get the angular momentum, and propagate a half step

      integrableObjects[i]->getJ( ji );

      for (j=0; j < 3; j++)
        ji[j] += dt2 * (Tb[j] * eConvert - ji[j]*chi);

      this->rotationPropagation( integrableObjects[i], ji );

      integrableObjects[i]->setJ( ji );
    }
  }

  // evolve chi and eta  half step

  evolveChiA();
  evolveEtaA();

  //calculate the integral of chidt
  integralOfChidt += dt2*chi;

  //save the old positions
  for(i = 0; i < integrableObjects.size(); i++){
    integrableObjects[i]->getPos(pos);
    for(j = 0; j < 3; j++)
      oldPos[i*3 + j] = pos[j];
  }

  //the first estimation of r(t+dt) is equal to  r(t)

  for(k = 0; k < 5; k ++){

    for(i =0 ; i < integrableObjects.size(); i++){

      integrableObjects[i]->getVel(vel);
      integrableObjects[i]->getPos(pos);

      this->getPosScale( pos, COM, i, sc );

      for(j = 0; j < 3; j++)
        pos[j] = oldPos[i*3 + j] + dt*(vel[j] + sc[j]);

      integrableObjects[i]->setPos( pos );
    }
    
    if(nConstrained)
      constrainA();
  }


  // Scale the box after all the positions have been moved:

  this->scaleSimBox();
}

template<typename T> void NPT<T>::moveB( void ){

  //new version of NPT
  int i, j, k;
  double Tb[3], ji[3], sc[3];
  double vel[3], frc[3];
  double mass;

  // Set things up for the iteration:

  for( i=0; i<integrableObjects.size(); i++ ){

    integrableObjects[i]->getVel( vel );

    for (j=0; j < 3; j++)
      oldVel[3*i + j]  = vel[j];

    if( integrableObjects[i]->isDirectional() ){

      integrableObjects[i]->getJ( ji );

      for (j=0; j < 3; j++)
        oldJi[3*i + j] = ji[j];

    }
  }

  // do the iteration:

  instaVol = tStats->getVolume();

  for (k=0; k < 4; k++) {

    instaTemp = tStats->getTemperature();
    instaPress = tStats->getPressure();

    // evolve chi another half step using the temperature at t + dt/2

    this->evolveChiB();
    this->evolveEtaB();
    this->calcVelScale();

    for( i=0; i<integrableObjects.size(); i++ ){

      integrableObjects[i]->getFrc( frc );
      integrableObjects[i]->getVel(vel);

      mass = integrableObjects[i]->getMass();

      getVelScaleB( sc, i );

      // velocity half step
      for (j=0; j < 3; j++)
        vel[j] = oldVel[3*i+j] + dt2 * ((frc[j] / mass ) * eConvert - sc[j]);

      integrableObjects[i]->setVel( vel );

      if( integrableObjects[i]->isDirectional() ){

        // get and convert the torque to body frame

        integrableObjects[i]->getTrq( Tb );
        integrableObjects[i]->lab2Body( Tb );

        for (j=0; j < 3; j++)
          ji[j] = oldJi[3*i + j] + dt2 * (Tb[j] * eConvert - oldJi[3*i+j]*chi);

          integrableObjects[i]->setJ( ji );
      }
    }

    if(nConstrained)
      constrainB();

    if ( this->chiConverged() && this->etaConverged() ) break;
  }

  //calculate integral of chida
  integralOfChidt += dt2*chi;


}

template<typename T> void NPT<T>::resetIntegrator() {
  chi = 0.0;
  T::resetIntegrator();
}

template<typename T> void NPT<T>::evolveChiA() {
  chi += dt2 * ( instaTemp / targetTemp - 1.0) / tt2;
  oldChi = chi;
}

template<typename T> void NPT<T>::evolveChiB() {

  prevChi = chi;
  chi = oldChi + dt2 * ( instaTemp / targetTemp - 1.0) / tt2;
}

template<typename T> bool NPT<T>::chiConverged() {

  return ( fabs( prevChi - chi ) <= chiTolerance );
}

template<typename T> int NPT<T>::readyCheck() {

  //check parent's readyCheck() first
  if (T::readyCheck() == -1)
    return -1;

  // First check to see if we have a target temperature.
  // Not having one is fatal.

  if (!have_target_temp) {
    sprintf( painCave.errMsg,
             "NPT error: You can't use the NPT integrator\n"
             "   without a targetTemp!\n"
             );
    painCave.isFatal = 1;
    simError();
    return -1;
  }

  if (!have_target_pressure) {
    sprintf( painCave.errMsg,
             "NPT error: You can't use the NPT integrator\n"
             "   without a targetPressure!\n"
             );
    painCave.isFatal = 1;
    simError();
    return -1;
  }

  // We must set tauThermostat.

  if (!have_tau_thermostat) {
    sprintf( painCave.errMsg,
             "NPT error: If you use the NPT\n"
             "   integrator, you must set tauThermostat.\n");
    painCave.isFatal = 1;
    simError();
    return -1;
  }

  // We must set tauBarostat.

  if (!have_tau_barostat) {
    sprintf( painCave.errMsg,
             "If you use the NPT integrator, you must set tauBarostat.\n");
    painCave.severity = OOPSE_ERROR;
    painCave.isFatal = 1;
    simError();
    return -1;
  }

  if (!have_chi_tolerance) {
    sprintf( painCave.errMsg,
             "Setting chi tolerance to 1e-6 in NPT integrator\n");
    chiTolerance = 1e-6;
    have_chi_tolerance = 1;
    painCave.severity = OOPSE_INFO;
    painCave.isFatal = 0;
    simError();
  }

  if (!have_eta_tolerance) {
    sprintf( painCave.errMsg,
             "Setting eta tolerance to 1e-6 in NPT integrator");
    etaTolerance = 1e-6;
    have_eta_tolerance = 1;
    painCave.severity = OOPSE_INFO;
    painCave.isFatal = 0;
    simError();
  }

  // We need NkBT a lot, so just set it here: This is the RAW number
  // of integrableObjects, so no subtraction or addition of constraints or
  // orientational degrees of freedom:

  NkBT = (double)(info->getTotIntegrableObjects()) * kB * targetTemp;

  // fkBT is used because the thermostat operates on more degrees of freedom
  // than the barostat (when there are particles with orientational degrees
  // of freedom).  

  fkBT = (double)(info->getNDF()) * kB * targetTemp;

  tt2 = tauThermostat * tauThermostat;
  tb2 = tauBarostat * tauBarostat;

  return 1;
}
