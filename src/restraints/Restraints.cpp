// Thermodynamic integration is not multiprocessor friendly right now

#include <iostream>
#include <stdlib.h>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <math.h>

using namespace std;

#include "restraints/Restraints.hpp"
#include "brains/SimInfo.hpp"
#include "utils/simError.h"
#include "io/basic_ifstrstream.hpp"

#ifdef IS_MPI
#include<mpi.h>
#include "brains/mpiSimulation.hpp"
#endif // is_mpi

#define PI 3.14159265359
#define TWO_PI 6.28318530718

Restraints::Restraints(double lambdaVal, double lambdaExp){
  lambdaValue = lambdaVal;
  lambdaK = lambdaExp;
  vector<double> resConsts;
  const char *jolt = " \t\n;,";

#ifdef IS_MPI
  if(worldRank == 0 ){
#endif // is_mpi

    strcpy(springName, "HarmSpringConsts.txt");
    
    ifstream springs(springName);
    
    if (!springs) { 
      sprintf(painCave.errMsg,
	      "Unable to open HarmSpringConsts.txt for reading, so the\n"
	      "\tdefault spring constants will be loaded. If you want\n"
	      "\tto specify spring constants, include a three line\n"
	      "\tHarmSpringConsts.txt file in the execution directory.\n");
      painCave.severity = OOPSE_WARNING;
      painCave.isFatal = 0;
      simError();   
      
      // load default spring constants
      kDist  = 6;  // spring constant in units of kcal/(mol*ang^2)
      kTheta = 7.5;   // in units of kcal/mol
      kOmega = 13.5;   // in units of kcal/mol
    } else  {
      
      springs.getline(inLine,999,'\n');
      // the file is blank!
      if (springs.eof()){
      sprintf(painCave.errMsg,
	      "HarmSpringConsts.txt file is not valid.\n"
	      "\tThe file should contain four rows, the last three containing\n"
	      "\ta label and the spring constant value. They should be listed\n"
	      "\tin the following order: kDist (positional restrant), kTheta\n"
	      "\t(rot. restraint: deflection of z-axis), and kOmega (rot.\n"
	      "\trestraint: rotation about the z-axis).\n");
      painCave.severity = OOPSE_ERROR;
      painCave.isFatal = 1;
      simError();   
      }
      // read in spring constants and check to make sure it is a valid file
      springs.getline(inLine,999,'\n');
      while (!springs.eof()){
	if (NULL != inLine){
	  token = strtok(inLine,jolt);
	  token = strtok(NULL,jolt);
	  if (NULL != token){
	    strcpy(inValue,token);
	    resConsts.push_back(atof(inValue));
	  }
 	}
	springs.getline(inLine,999,'\n');
      }
      if (resConsts.size() == 3){
	kDist = resConsts[0];
	kTheta = resConsts[1];
	kOmega = resConsts[2];
      }
      else {
	sprintf(painCave.errMsg,
		"HarmSpringConsts.txt file is not valid.\n"
		"\tThe file should contain four rows, the last three containing\n"
		"\ta label and the spring constant value. They should be listed\n"
		"\tin the following order: kDist (positional restrant), kTheta\n"
		"\t(rot. restraint: deflection of z-axis), and kOmega (rot.\n"
		"\trestraint: rotation about the z-axis).\n");
	painCave.severity = OOPSE_ERROR;
	painCave.isFatal = 1;
	simError();  
      }
    }
#ifdef IS_MPI
  }
  
  MPI_Bcast(&kDist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&kTheta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&kOmega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
  
  sprintf( checkPointMsg,
           "Sucessfully opened and read spring file.\n");
  MPIcheckPoint();

#endif // is_mpi
  
  sprintf(painCave.errMsg,
	  "The spring constants for thermodynamic integration are:\n"
	  "\tkDist = %lf\n"
	  "\tkTheta = %lf\n"
	  "\tkOmega = %lf\n", kDist, kTheta, kOmega);
  painCave.severity = OOPSE_INFO;
  painCave.isFatal = 0;
  simError();   
}

Restraints::~Restraints(){
}

void Restraints::Calc_rVal(double position[3], double refPosition[3]){
  delRx = position[0] - refPosition[0];
  delRy = position[1] - refPosition[1];
  delRz = position[2] - refPosition[2];

  return;
}

void Restraints::Calc_body_thetaVal(double matrix[3][3], double refUnit[3]){
  ub0x = matrix[0][0]*refUnit[0] + matrix[0][1]*refUnit[1]
    + matrix[0][2]*refUnit[2];
  ub0y = matrix[1][0]*refUnit[0] + matrix[1][1]*refUnit[1]
    + matrix[1][2]*refUnit[2];
  ub0z = matrix[2][0]*refUnit[0] + matrix[2][1]*refUnit[1]
    + matrix[2][2]*refUnit[2];

  normalize = sqrt(ub0x*ub0x + ub0y*ub0y + ub0z*ub0z);
  ub0x = ub0x/normalize;
  ub0y = ub0y/normalize;
  ub0z = ub0z/normalize;

  // Theta is the dot product of the reference and new z-axes
  theta = acos(ub0z);

  return;
}

void Restraints::Calc_body_omegaVal(double matrix[3][3], double zAngle){
  double zRotator[3][3];
  double tempOmega;
  double wholeTwoPis;
  // Use the omega accumulated from the rotation propagation
  omega = zAngle;

  // translate the omega into a range between -PI and PI
  if (omega < -PI){
    tempOmega = omega / -TWO_PI;
    wholeTwoPis = floor(tempOmega);
    tempOmega = omega + TWO_PI*wholeTwoPis;
    if (tempOmega < -PI)
      omega = tempOmega + TWO_PI;
    else
      omega = tempOmega;
  }
  if (omega > PI){
    tempOmega = omega / TWO_PI;
    wholeTwoPis = floor(tempOmega);
    tempOmega = omega - TWO_PI*wholeTwoPis;
    if (tempOmega > PI)
      omega = tempOmega - TWO_PI;   
    else
      omega = tempOmega;
  }

  vb0x = sin(omega);
  vb0y = cos(omega);
  vb0z = 0.0;

  normalize = sqrt(vb0x*vb0x + vb0y*vb0y + vb0z*vb0z);
  vb0x = vb0x/normalize;
  vb0y = vb0y/normalize;
  vb0z = vb0z/normalize;

  return;
}

double Restraints::Calc_Restraint_Forces(vector<StuntDouble*> vecParticles){
  double pos[3];
  double A[3][3];
  double refPos[3];
  double refVec[3];
  double tolerance;
  double tempPotent;
  double factor;
  double spaceTrq[3];
  double omegaPass;
  GenericData* data;
  DoubleGenericData* doubleData;

  tolerance = 5.72957795131e-7;

  harmPotent = 0.0;  // zero out the global harmonic potential variable

  factor = 1 - pow(lambdaValue, lambdaK);

  for (i=0; i<vecParticles.size(); i++){
    // obtain the current and reference positions
    vecParticles[i]->getPos(pos);

    data = vecParticles[i]->getProperty("refPosX");
    if (data){
      doubleData = dynamic_cast<DoubleGenericData*>(data);
      if (!doubleData){
	cerr << "Can't obtain refPosX from StuntDouble\n";
	return 0.0;
      }
      else refPos[0] = doubleData->getData();
    }
    data = vecParticles[i]->getProperty("refPosY");
    if (data){
      doubleData = dynamic_cast<DoubleGenericData*>(data);
      if (!doubleData){
	cerr << "Can't obtain refPosY from StuntDouble\n";
	return 0.0;
      }
      else refPos[1] = doubleData->getData();
    }
    data = vecParticles[i]->getProperty("refPosZ");
    if (data){
      doubleData = dynamic_cast<DoubleGenericData*>(data);
      if (!doubleData){
	cerr << "Can't obtain refPosZ from StuntDouble\n";
	return 0.0;
      }
      else refPos[2] = doubleData->getData();
    }

    // calculate the displacement
    Calc_rVal( pos, refPos );

    // calculate the derivatives
    dVdrx = -kDist*delRx;
    dVdry = -kDist*delRy;
    dVdrz = -kDist*delRz;
    
    // next we calculate the restraint forces
    restraintFrc[0] = dVdrx;
    restraintFrc[1] = dVdry;
    restraintFrc[2] = dVdrz;
    tempPotent = 0.5*kDist*(delRx*delRx + delRy*delRy + delRz*delRz);

    // apply the lambda scaling factor to the forces
    for (j = 0; j < 3; j++) restraintFrc[j] *= factor;

    // and add the temporary force to the total force
    vecParticles[i]->addFrc(restraintFrc);

    // if the particle is directional, we accumulate the rot. restraints
    if (vecParticles[i]->isDirectional()){
 
      // get the current rotation matrix and reference vector
      vecParticles[i]->getA(A);
      
      data = vecParticles[i]->getProperty("refVectorX");
      if (data){
	doubleData = dynamic_cast<DoubleGenericData*>(data);
	if (!doubleData){
	  cerr << "Can't obtain refVectorX from StuntDouble\n";
	  return 0.0;
	}
	else refVec[0] = doubleData->getData();
      }
      data = vecParticles[i]->getProperty("refVectorY");
      if (data){
	doubleData = dynamic_cast<DoubleGenericData*>(data);
	if (!doubleData){
	  cerr << "Can't obtain refVectorY from StuntDouble\n";
	  return 0.0;
	}
	else refVec[1] = doubleData->getData();
      }
      data = vecParticles[i]->getProperty("refVectorZ");
      if (data){
	doubleData = dynamic_cast<DoubleGenericData*>(data);
	if (!doubleData){
	  cerr << "Can't obtain refVectorZ from StuntDouble\n";
	  return 0.0;
	}
	else refVec[2] = doubleData->getData();
      }
      
      // calculate the theta and omega displacements
      Calc_body_thetaVal( A, refVec );
      omegaPass = vecParticles[i]->getZangle();
      Calc_body_omegaVal( A, omegaPass );

      // uTx... and vTx... are the body-fixed z and y unit vectors
      uTx = 0.0;
      uTy = 0.0;
      uTz = 1.0;
      vTx = 0.0;
      vTy = 1.0;
      vTz = 0.0;

      dVdux = 0.0;
      dVduy = 0.0;
      dVduz = 0.0;
      dVdvx = 0.0;
      dVdvy = 0.0;
      dVdvz = 0.0;

      if (fabs(theta) > tolerance) {
	dVdux = -(kTheta*theta/sin(theta))*ub0x;
	dVduy = -(kTheta*theta/sin(theta))*ub0y;
	dVduz = -(kTheta*theta/sin(theta))*ub0z;
      }

      if (fabs(omega) > tolerance) {
	dVdvx = -(kOmega*omega/sin(omega))*vb0x;
	dVdvy = -(kOmega*omega/sin(omega))*vb0y;
	dVdvz = -(kOmega*omega/sin(omega))*vb0z;
      }

      // next we calculate the restraint torques
      restraintTrq[0] = 0.0;
      restraintTrq[1] = 0.0;
      restraintTrq[2] = 0.0;

      if (fabs(omega) > tolerance) {
	restraintTrq[0] += 0.0;
	restraintTrq[1] += 0.0;
	restraintTrq[2] += vTy*dVdvx;
	tempPotent += 0.5*(kOmega*omega*omega);
      }
      if (fabs(theta) > tolerance) {
	restraintTrq[0] += (uTz*dVduy);
	restraintTrq[1] += -(uTz*dVdux);
	restraintTrq[2] += 0.0;
	tempPotent += 0.5*(kTheta*theta*theta);
      }

      // apply the lambda scaling factor to these torques
      for (j = 0; j < 3; j++) restraintTrq[j] *= factor;

      // now we need to convert from body-fixed torques to space-fixed torques
      spaceTrq[0] = A[0][0]*restraintTrq[0] + A[1][0]*restraintTrq[1] 
	+ A[2][0]*restraintTrq[2];
      spaceTrq[1] = A[0][1]*restraintTrq[0] + A[1][1]*restraintTrq[1] 
	+ A[2][1]*restraintTrq[2];
      spaceTrq[2] = A[0][2]*restraintTrq[0] + A[1][2]*restraintTrq[1] 
	+ A[2][2]*restraintTrq[2];

      // now pass this temporary torque vector to the total torque
      vecParticles[i]->addTrq(spaceTrq);
    }

    // update the total harmonic potential with this object's contribution
    harmPotent += tempPotent;
  }
  
  // we can finish by returning the appropriately scaled potential energy
  tempPotent = harmPotent * factor;
  return tempPotent;
}

void Restraints::Write_zAngle_File(vector<StuntDouble*> vecParticles,
				   int currTime,
				   int nIntObj){

  char zOutName[200];

  std::cerr << nIntObj << " is the number of integrable objects\n";

  //#ifndef IS_MPI
  
  strcpy(zOutName,"zAngle.ang");
  
  ofstream angleOut(zOutName);
  angleOut << currTime << ": omega values at this time\n";
  for (i=0; i<vecParticles.size(); i++) {
    angleOut << vecParticles[i]->getZangle() << "\n";
  }

  return;
}

double Restraints::getVharm(){
  return harmPotent;
}

