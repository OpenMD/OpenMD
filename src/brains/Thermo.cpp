#include <math.h>
#include <iostream>
using namespace std;

#ifdef IS_MPI
#include <mpi.h>
#endif //is_mpi

#include "Thermo.hpp"
#include "SRI.hpp"
#include "Integrator.hpp"
#include "simError.h"
#include "MatVec3.h"

#ifdef IS_MPI
#define __C
#include "mpiSimulation.hpp"
#endif // is_mpi

inline double roundMe( double x ){
	  return ( x >= 0 ) ? floor( x + 0.5 ) : ceil( x - 0.5 );
}

Thermo::Thermo( SimInfo* the_info ) { 
  info = the_info;
  int baseSeed = the_info->getSeed();
  
  gaussStream = new gaussianSPRNG( baseSeed );
}

Thermo::~Thermo(){
  delete gaussStream;
}

double Thermo::getKinetic(){

  const double e_convert = 4.184E-4; // convert kcal/mol -> (amu A^2)/fs^2
  double kinetic;
  double amass;
  double aVel[3], aJ[3], I[3][3];
  int i, j, k, kl;

  double kinetic_global;
  vector<StuntDouble *> integrableObjects = info->integrableObjects;
  
  kinetic = 0.0;
  kinetic_global = 0.0;

  for (kl=0; kl<integrableObjects.size(); kl++) {
    integrableObjects[kl]->getVel(aVel);
    amass = integrableObjects[kl]->getMass();

   for(j=0; j<3; j++) 
      kinetic += amass*aVel[j]*aVel[j];

   if (integrableObjects[kl]->isDirectional()){
 
      integrableObjects[kl]->getJ( aJ );
      integrableObjects[kl]->getI( I );

      if (integrableObjects[kl]->isLinear()) {
        i = integrableObjects[kl]->linearAxis();
        j = (i+1)%3;
        k = (i+2)%3;
        kinetic += aJ[j]*aJ[j]/I[j][j] + aJ[k]*aJ[k]/I[k][k];
      } else {
        for (j=0; j<3; j++) 
          kinetic += aJ[j]*aJ[j] / I[j][j];
      }
   }
  }
#ifdef IS_MPI
  MPI_Allreduce(&kinetic,&kinetic_global,1,MPI_DOUBLE,
		MPI_SUM, MPI_COMM_WORLD);
  kinetic = kinetic_global;
#endif //is_mpi
  
  kinetic = kinetic * 0.5 / e_convert;

  return kinetic;
}

double Thermo::getPotential(){
  
  double potential_local;
  double potential;
  int el, nSRI;
  Molecule* molecules;

  molecules = info->molecules;
  nSRI = info->n_SRI;

  potential_local = 0.0;
  potential = 0.0;
  potential_local += info->lrPot;

  for( el=0; el<info->n_mol; el++ ){    
    potential_local += molecules[el].getPotential();
  }

  // Get total potential for entire system from MPI.
#ifdef IS_MPI
  MPI_Allreduce(&potential_local,&potential,1,MPI_DOUBLE,
		MPI_SUM, MPI_COMM_WORLD);
#else
  potential = potential_local; 
#endif // is_mpi

  return potential;
}

double Thermo::getTotalE(){

  double total;

  total = this->getKinetic() + this->getPotential();
  return total;
}

double Thermo::getTemperature(){

  const double kb = 1.9872156E-3; // boltzman's constant in kcal/(mol K)
  double temperature;

  temperature = ( 2.0 * this->getKinetic() ) / ((double)info->ndf * kb );
  return temperature;
}

double Thermo::getVolume() {

  return info->boxVol;
}

double Thermo::getPressure() {

  // Relies on the calculation of the full molecular pressure tensor
  
  const double p_convert = 1.63882576e8;
  double press[3][3];
  double pressure;

  this->getPressureTensor(press);

  pressure = p_convert * (press[0][0] + press[1][1] + press[2][2]) / 3.0;

  return pressure;
}

double Thermo::getPressureX() {

  // Relies on the calculation of the full molecular pressure tensor
  
  const double p_convert = 1.63882576e8;
  double press[3][3];
  double pressureX;

  this->getPressureTensor(press);

  pressureX = p_convert * press[0][0];

  return pressureX;
}

double Thermo::getPressureY() {

  // Relies on the calculation of the full molecular pressure tensor
  
  const double p_convert = 1.63882576e8;
  double press[3][3];
  double pressureY;

  this->getPressureTensor(press);

  pressureY = p_convert * press[1][1];

  return pressureY;
}

double Thermo::getPressureZ() {

  // Relies on the calculation of the full molecular pressure tensor
  
  const double p_convert = 1.63882576e8;
  double press[3][3];
  double pressureZ;

  this->getPressureTensor(press);

  pressureZ = p_convert * press[2][2];

  return pressureZ;
}


void Thermo::getPressureTensor(double press[3][3]){
  // returns pressure tensor in units amu*fs^-2*Ang^-1
  // routine derived via viral theorem description in:
  // Paci, E. and Marchi, M. J.Phys.Chem. 1996, 100, 4314-4322

  const double e_convert = 4.184e-4;

  double molmass, volume;
  double vcom[3];
  double p_local[9], p_global[9];
  int i, j, k;

  for (i=0; i < 9; i++) {    
    p_local[i] = 0.0;
    p_global[i] = 0.0;
  }

  // use velocities of integrableObjects and their masses:  

  for (i=0; i < info->integrableObjects.size(); i++) {

    molmass = info->integrableObjects[i]->getMass();
    
    info->integrableObjects[i]->getVel(vcom);
    
    p_local[0] += molmass * (vcom[0] * vcom[0]); 
    p_local[1] += molmass * (vcom[0] * vcom[1]); 
    p_local[2] += molmass * (vcom[0] * vcom[2]); 
    p_local[3] += molmass * (vcom[1] * vcom[0]); 
    p_local[4] += molmass * (vcom[1] * vcom[1]); 
    p_local[5] += molmass * (vcom[1] * vcom[2]); 
    p_local[6] += molmass * (vcom[2] * vcom[0]); 
    p_local[7] += molmass * (vcom[2] * vcom[1]); 
    p_local[8] += molmass * (vcom[2] * vcom[2]); 

  }

  // Get total for entire system from MPI.
  
#ifdef IS_MPI
  MPI_Allreduce(p_local,p_global,9,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  for (i=0; i<9; i++) {
    p_global[i] = p_local[i]; 
  }
#endif // is_mpi

  volume = this->getVolume();



  for(i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      k = 3*i + j;
      press[i][j] = (p_global[k] + info->tau[k]*e_convert) / volume;
    }
  }
}

void Thermo::velocitize() {
  
  double aVel[3], aJ[3], I[3][3];
  int i, j, l, m, n, vr, vd; // velocity randomizer loop counters
  double vdrift[3];
  double vbar;
  const double kb = 8.31451e-7; // kb in amu, angstroms, fs, etc.
  double av2;
  double kebar;
  double temperature;
  int nobj;

  if (!info->have_target_temp) {
    sprintf( painCave.errMsg,
             "You can't resample the velocities without a targetTemp!\n"
             );
    painCave.isFatal = 1;
    painCave.severity = OOPSE_ERROR;
    simError();
    return;
  }

  nobj = info->integrableObjects.size();
  
  temperature   = info->target_temp;
  
  kebar = kb * temperature * (double)info->ndfRaw / 
    ( 2.0 * (double)info->ndf );
  
  for(vr = 0; vr < nobj; vr++){
    
    // uses equipartition theory to solve for vbar in angstrom/fs

    av2 = 2.0 * kebar / info->integrableObjects[vr]->getMass();
    vbar = sqrt( av2 );

    // picks random velocities from a gaussian distribution
    // centered on vbar

    for (j=0; j<3; j++) 
      aVel[j] = vbar * gaussStream->getGaussian();
    
    info->integrableObjects[vr]->setVel( aVel );
    
    if(info->integrableObjects[vr]->isDirectional()){

      info->integrableObjects[vr]->getI( I );

      if (info->integrableObjects[vr]->isLinear()) {

        l= info->integrableObjects[vr]->linearAxis();
        m = (l+1)%3;
        n = (l+2)%3;

        aJ[l] = 0.0;
        vbar = sqrt( 2.0 * kebar * I[m][m] );
        aJ[m] = vbar * gaussStream->getGaussian();
        vbar = sqrt( 2.0 * kebar * I[n][n] );
        aJ[n] = vbar * gaussStream->getGaussian();
        
      } else {
        for (j = 0 ; j < 3; j++) {
          vbar = sqrt( 2.0 * kebar * I[j][j] );
          aJ[j] = vbar * gaussStream->getGaussian();
        }	
      } // else isLinear

      info->integrableObjects[vr]->setJ( aJ ); 
      
    }//isDirectional 

  }

  // Get the Center of Mass drift velocity.

  getCOMVel(vdrift);
  
  //  Corrects for the center of mass drift.
  // sums all the momentum and divides by total mass.

  for(vd = 0; vd < nobj; vd++){
    
    info->integrableObjects[vd]->getVel(aVel);
    
    for (j=0; j < 3; j++) 
      aVel[j] -= vdrift[j];
        
    info->integrableObjects[vd]->setVel( aVel );
  }

}

void Thermo::getCOMVel(double vdrift[3]){

  double mtot, mtot_local;
  double aVel[3], amass;
  double vdrift_local[3];
  int vd, j;
  int nobj;

  nobj   = info->integrableObjects.size();

  mtot_local = 0.0;
  vdrift_local[0] = 0.0;
  vdrift_local[1] = 0.0;
  vdrift_local[2] = 0.0;
  
  for(vd = 0; vd < nobj; vd++){
    
    amass = info->integrableObjects[vd]->getMass();
    info->integrableObjects[vd]->getVel( aVel );

    for(j = 0; j < 3; j++) 
      vdrift_local[j] += aVel[j] * amass;
    
    mtot_local += amass;
  }

#ifdef IS_MPI
  MPI_Allreduce(&mtot_local,&mtot,1,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(vdrift_local,vdrift,3,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
#else
  mtot = mtot_local;
  for(vd = 0; vd < 3; vd++) {
    vdrift[vd] = vdrift_local[vd];
  }
#endif
    
  for (vd = 0; vd < 3; vd++) {
    vdrift[vd] = vdrift[vd] / mtot;
  }
  
}

void Thermo::getCOM(double COM[3]){

  double mtot, mtot_local;
  double aPos[3], amass;
  double COM_local[3];
  int i, j;
  int nobj;

  mtot_local = 0.0;
  COM_local[0] = 0.0;
  COM_local[1] = 0.0;
  COM_local[2] = 0.0;

  nobj = info->integrableObjects.size();
  for(i = 0; i < nobj; i++){
    
    amass = info->integrableObjects[i]->getMass();
    info->integrableObjects[i]->getPos( aPos );

    for(j = 0; j < 3; j++) 
      COM_local[j] += aPos[j] * amass;
    
    mtot_local += amass;
  }

#ifdef IS_MPI
  MPI_Allreduce(&mtot_local,&mtot,1,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(COM_local,COM,3,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
#else
  mtot = mtot_local;
  for(i = 0; i < 3; i++) {
    COM[i] = COM_local[i];
  }
#endif
    
  for (i = 0; i < 3; i++) {
    COM[i] = COM[i] / mtot;
  }
}

void Thermo::removeCOMdrift() {
  double vdrift[3], aVel[3];
  int vd, j, nobj;

  nobj = info->integrableObjects.size();

  // Get the Center of Mass drift velocity.

  getCOMVel(vdrift);
  
  //  Corrects for the center of mass drift.
  // sums all the momentum and divides by total mass.

  for(vd = 0; vd < nobj; vd++){
    
    info->integrableObjects[vd]->getVel(aVel);
    
    for (j=0; j < 3; j++) 
      aVel[j] -= vdrift[j];
        
    info->integrableObjects[vd]->setVel( aVel );
  }
}
