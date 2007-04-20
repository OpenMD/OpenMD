/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 */

#include <stdlib.h>
#include <math.h>

using namespace std;

#include "restraints/Restraints.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"

#define PI 3.14159265359
#define TWO_PI 6.28318530718

namespace oopse {
  
  Restraints::Restraints(SimInfo* info, RealType lambdaVal, RealType lambdaExp){
    info_ = info;
    Globals* simParam = info_->getSimParams();

    lambdaValue = lambdaVal;
    lambdaK = lambdaExp;
    
    if (simParam->getUseSolidThermInt()) {
      if (simParam->haveThermIntDistSpringConst()) {
        kDist = simParam->getThermIntDistSpringConst();
      }
      else{
        kDist = simParam->getThermIntDistSpringConst();
        sprintf(painCave.errMsg,
                "ThermoIntegration Warning: the spring constant for the\n"
                "\ttranslational restraint was not specified. OOPSE will use\n"
                "\ta default value of %f. To set it to something else, use\n"
                "\tthe thermIntDistSpringConst variable.\n",
                kDist);
        painCave.isFatal = 0;
        simError(); 
      }
      if (simParam->haveThermIntThetaSpringConst()) {
        kTheta = simParam->getThermIntThetaSpringConst();
      }
      else{
        kTheta = simParam->getThermIntThetaSpringConst();
        sprintf(painCave.errMsg,
                "ThermoIntegration Warning: the spring constant for the\n"
                "\tdeflection orientational restraint was not specified.\n"
                "\tOOPSE will use a default value of %f. To set it to\n"
                "\tsomething else, use the thermIntThetaSpringConst variable.\n",
                kTheta);
        painCave.isFatal = 0;
        simError(); 
      }
      if (simParam->haveThermIntOmegaSpringConst()) {
        kOmega = simParam->getThermIntOmegaSpringConst();
      }
      else{
        kOmega = simParam->getThermIntOmegaSpringConst();
        sprintf(painCave.errMsg,
                "ThermoIntegration Warning: the spring constant for the\n"
                "\tspin orientational restraint was not specified. OOPSE\n"
                "\twill use a default value of %f. To set it to something\n"
                "\telse, use the thermIntOmegaSpringConst variable.\n",
                kOmega);
        painCave.isFatal = 0;
        simError(); 
      }
    }
    
    // build a RestReader and read in important information
    
    restRead_ = new RestReader(info_);
    restRead_->readIdealCrystal();
    restRead_->readZangle();
    
    delete restRead_;
    restRead_ = NULL;
    
  }
  
  Restraints::~Restraints(){
  }
  
  void Restraints::Calc_rVal(Vector3d &position, RealType refPosition[3]){
    delRx = position.x() - refPosition[0];
    delRy = position.y() - refPosition[1];
    delRz = position.z() - refPosition[2];
    
    return;
  }
  
  void Restraints::Calc_body_thetaVal(RotMat3x3d &matrix, RealType refUnit[3]){
    ub0x = matrix(0,0)*refUnit[0] + matrix(0,1)*refUnit[1]
      + matrix(0,2)*refUnit[2];
    ub0y = matrix(1,0)*refUnit[0] + matrix(1,1)*refUnit[1]
      + matrix(1,2)*refUnit[2];
    ub0z = matrix(2,0)*refUnit[0] + matrix(2,1)*refUnit[1]
      + matrix(2,2)*refUnit[2];
    
    normalize = sqrt(ub0x*ub0x + ub0y*ub0y + ub0z*ub0z);
    ub0x = ub0x/normalize;
    ub0y = ub0y/normalize;
    ub0z = ub0z/normalize;
    
    // Theta is the dot product of the reference and new z-axes
    theta = acos(ub0z);
    
    return;
  }
  
  void Restraints::Calc_body_omegaVal(RealType zAngle){
    RealType tempOmega;
    RealType wholeTwoPis;
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
  
  RealType Restraints::Calc_Restraint_Forces(){
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::IntegrableObjectIterator ii;
    StuntDouble* integrableObject;
    Vector3d pos;
    RotMat3x3d A;
    RealType refPos[3];
    RealType refVec[3];
    RealType tolerance;
    RealType tempPotent;
    RealType factor;
    RealType spaceTrq[3];
    RealType omegaPass;
    GenericData* data;
    DoubleGenericData* doubleData;
    
    tolerance = 5.72957795131e-7;
    
    harmPotent = 0.0;  // zero out the global harmonic potential variable
    
    factor = 1 - pow(lambdaValue, lambdaK);
    
    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) {
      for (integrableObject = mol->beginIntegrableObject(ii); 
           integrableObject != NULL; 
           integrableObject = mol->nextIntegrableObject(ii)) {
        
        // obtain the current and reference positions
        pos = integrableObject->getPos();
        
        data = integrableObject->getPropertyByName("refPosX");
        if (data){
          doubleData = dynamic_cast<DoubleGenericData*>(data);
          if (!doubleData){
            cerr << "Can't obtain refPosX from StuntDouble\n";
            return 0.0;
          }
          else refPos[0] = doubleData->getData();
        }
        data = integrableObject->getPropertyByName("refPosY");
        if (data){
          doubleData = dynamic_cast<DoubleGenericData*>(data);
          if (!doubleData){
            cerr << "Can't obtain refPosY from StuntDouble\n";
            return 0.0;
          }
          else refPos[1] = doubleData->getData();
        }
        data = integrableObject->getPropertyByName("refPosZ");
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
        integrableObject->addFrc(restraintFrc);
        
        // if the particle is directional, we accumulate the rot. restraints
        if (integrableObject->isDirectional()){
          
          // get the current rotation matrix and reference vector
          A = integrableObject->getA();
          
          data = integrableObject->getPropertyByName("refVectorX");
          if (data){
            doubleData = dynamic_cast<DoubleGenericData*>(data);
            if (!doubleData){
              cerr << "Can't obtain refVectorX from StuntDouble\n";
              return 0.0;
            }
            else refVec[0] = doubleData->getData();
          }
          data = integrableObject->getPropertyByName("refVectorY");
          if (data){
            doubleData = dynamic_cast<DoubleGenericData*>(data);
            if (!doubleData){
              cerr << "Can't obtain refVectorY from StuntDouble\n";
              return 0.0;
            }
            else refVec[1] = doubleData->getData();
          }
          data = integrableObject->getPropertyByName("refVectorZ");
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
          omegaPass = integrableObject->getZangle();
          Calc_body_omegaVal( omegaPass );
          
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
          
          // now we need to convert from body-fixed to space-fixed torques
          spaceTrq[0] = A(0,0)*restraintTrq[0] + A(1,0)*restraintTrq[1] 
            + A(2,0)*restraintTrq[2];
          spaceTrq[1] = A(0,1)*restraintTrq[0] + A(1,1)*restraintTrq[1] 
            + A(2,1)*restraintTrq[2];
          spaceTrq[2] = A(0,2)*restraintTrq[0] + A(1,2)*restraintTrq[1] 
            + A(2,2)*restraintTrq[2];
          
          // now pass this temporary torque vector to the total torque
          integrableObject->addTrq(spaceTrq);
        }
        
        // update the total harmonic potential with this object's contribution
        harmPotent += tempPotent;
      }
      
    }
    
    // we can finish by returning the appropriately scaled potential energy
    tempPotent = harmPotent * factor;
    
    return tempPotent;
    
  }
  
}// end namespace oopse
