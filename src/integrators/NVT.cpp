#include <math.h>

#include "Atom.hpp"
#include "SRI.hpp"
#include "AbstractClasses.hpp"
#include "SimInfo.hpp"
#include "ForceFields.hpp"
#include "Thermo.hpp"
#include "ReadWrite.hpp"
#include "Integrator.hpp"
#include "simError.h"


// Basic thermostating via Hoover, Phys.Rev.A, 1985, Vol. 31 (5) 1695-1697

template<typename T> NVT<T>::NVT ( SimInfo *theInfo, ForceFields* the_ff):
  T( theInfo, the_ff )
{
  GenericData* data;
  DoubleData * chiValue;
  DoubleData * integralOfChidtValue;

  chiValue = NULL;
  integralOfChidtValue = NULL;

  chi = 0.0;
  have_tau_thermostat = 0;
  have_target_temp = 0;
  have_chi_tolerance = 0;
  integralOfChidt = 0.0;


  if( theInfo->useInitXSstate ){

    // retrieve chi and integralOfChidt from simInfo
    data = info->getProperty(CHIVALUE_ID);
    if(data){
      chiValue = dynamic_cast<DoubleData*>(data);
    }
    
    data = info->getProperty(INTEGRALOFCHIDT_ID);
    if(data){
      integralOfChidtValue = dynamic_cast<DoubleData*>(data);
    }
    
    // chi and integralOfChidt should appear by pair
    if(chiValue && integralOfChidtValue){
      chi = chiValue->getData();
      integralOfChidt = integralOfChidtValue->getData();
    }
  }

  oldVel = new double[3*integrableObjects.size()];
  oldJi = new double[3*integrableObjects.size()];
}

template<typename T> NVT<T>::~NVT() {
  delete[] oldVel;
  delete[] oldJi;
}

template<typename T> void NVT<T>::moveA() {

  int i, j;
  DirectionalAtom* dAtom;
  double Tb[3], ji[3];
  double mass;
  double vel[3], pos[3], frc[3];

  double instTemp;

  // We need the temperature at time = t for the chi update below:

  instTemp = tStats->getTemperature();

  for( i=0; i < integrableObjects.size(); i++ ){

    integrableObjects[i]->getVel( vel );
    integrableObjects[i]->getPos( pos );
    integrableObjects[i]->getFrc( frc );

    mass = integrableObjects[i]->getMass();

    for (j=0; j < 3; j++) {
      // velocity half step  (use chi from previous step here):
      vel[j] += dt2 * ((frc[j] / mass ) * eConvert - vel[j]*chi);
      // position whole step
      pos[j] += dt * vel[j];
    }

    integrableObjects[i]->setVel( vel );
    integrableObjects[i]->setPos( pos );

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
  
  if(nConstrained)
    constrainA();

  // Finally, evolve chi a half step (just like a velocity) using
  // temperature at time t, not time t+dt/2

  //std::cerr << "targetTemp = " << targetTemp << " instTemp = " << instTemp << " tauThermostat = " << tauThermostat << " integral of Chi = " << integralOfChidt << "\n";
  
  chi += dt2 * ( instTemp / targetTemp - 1.0) / (tauThermostat*tauThermostat);
  integralOfChidt += chi*dt2;

}

template<typename T> void NVT<T>::moveB( void ){
  int i, j, k;
  double Tb[3], ji[3];
  double vel[3], frc[3];
  double mass;
  double instTemp;
  double oldChi, prevChi;

  // Set things up for the iteration:

  oldChi = chi;

  for( i=0; i < integrableObjects.size(); i++ ){

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

  for (k=0; k < 4; k++) {

    instTemp = tStats->getTemperature();

    // evolve chi another half step using the temperature at t + dt/2

    prevChi = chi;
    chi = oldChi + dt2 * ( instTemp / targetTemp - 1.0) /
      (tauThermostat*tauThermostat);

    for( i=0; i < integrableObjects.size(); i++ ){

      integrableObjects[i]->getFrc( frc );
      integrableObjects[i]->getVel(vel);

      mass = integrableObjects[i]->getMass();

      // velocity half step
      for (j=0; j < 3; j++)
	vel[j] = oldVel[3*i+j] + dt2 * ((frc[j] / mass ) * eConvert - oldVel[3*i + j]*chi);

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

    if (fabs(prevChi - chi) <= chiTolerance) break;
  }

  integralOfChidt += dt2*chi;
}

template<typename T> void NVT<T>::resetIntegrator( void ){

  chi = 0.0;
  integralOfChidt = 0.0;
}

template<typename T> int NVT<T>::readyCheck() {

  //check parent's readyCheck() first
  if (T::readyCheck() == -1)
    return -1;

  // First check to see if we have a target temperature.
  // Not having one is fatal.

  if (!have_target_temp) {
    sprintf( painCave.errMsg,
             "You can't use the NVT integrator without a targetTemp!\n"
             );
    painCave.isFatal = 1;
    painCave.severity = OOPSE_ERROR;
    simError();
    return -1;
  }

  // We must set tauThermostat.

  if (!have_tau_thermostat) {
    sprintf( painCave.errMsg,
             "If you use the constant temperature\n"
             "\tintegrator, you must set tauThermostat.\n");
    painCave.severity = OOPSE_ERROR;
    painCave.isFatal = 1;
    simError();
    return -1;
  }

  if (!have_chi_tolerance) {
    sprintf( painCave.errMsg,
             "In NVT integrator: setting chi tolerance to 1e-6\n");
    chiTolerance = 1e-6;
    have_chi_tolerance = 1;
    painCave.severity = OOPSE_INFO;
    painCave.isFatal = 0;
    simError();
  }

  return 1;

}

template<typename T> double NVT<T>::getConservedQuantity(void){

  double conservedQuantity;
  double fkBT;
  double Energy;
  double thermostat_kinetic;
  double thermostat_potential;

  fkBT = (double)(info->ndf) * kB * targetTemp;

  Energy = tStats->getTotalE();

  thermostat_kinetic = fkBT* tauThermostat * tauThermostat * chi * chi /
    (2.0 * eConvert);

  thermostat_potential = fkBT * integralOfChidt / eConvert;

  conservedQuantity = Energy + thermostat_kinetic + thermostat_potential;

  return conservedQuantity;
}

template<typename T> string NVT<T>::getAdditionalParameters(void){
  string parameters;
  const int BUFFERSIZE = 2000; // size of the read buffer
  char buffer[BUFFERSIZE];

  sprintf(buffer,"\t%G\t%G;", chi, integralOfChidt);
  parameters += buffer;

  return parameters;
}
