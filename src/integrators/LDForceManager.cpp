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
#include <fstream> 
#include "integrators/LDForceManager.hpp"
#include "math/CholeskyDecomposition.hpp"
#include "utils/OOPSEConstant.hpp"
namespace oopse {

  LDForceManager::LDForceManager(SimInfo* info) : ForceManager(info){
    Globals* simParams = info->getSimParams();
    std::map<std::string, HydroProp> hydroPropMap;
    if (simParams->haveHydroPropFile()) {
        hydroPropMap = parseFrictionFile(simParams->getHydroPropFile());
    } else {
        //error
    }

    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* integrableObject;
    for (mol = info->beginMolecule(i); mol != NULL; mol = info->nextMolecule(i)) {
      for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL;
	      integrableObject = mol->nextIntegrableObject(j)) {
            std::map<std::string, HydroProp>::iterator iter = hydroPropMap.find(integrableObject->getType());
            if (iter != hydroPropMap.end()) {
                hydroProps_.push_back(iter->second);
            } else {
                //error
            }
            
	   }
    }
    variance_ = 2.0 * OOPSEConstant::kb*simParams->getTargetTemp()/simParams->getDt();
  }
  std::map<std::string, HydroProp> LDForceManager::parseFrictionFile(const std::string& filename) {
    std::map<std::string, HydroProp> props;
    std::ifstream ifs(filename.c_str());
    if (ifs.is_open()) {

    }

    const unsigned int BufferSize = 65535;
    char buffer[BufferSize];   
    while (ifs.getline(buffer, BufferSize)) {
        StringTokenizer tokenizer(buffer);
        HydroProp currProp;
        if (tokenizer.countTokens() >= 67) {
            std::string atomName = tokenizer.nextToken();
            currProp.cod[0] = tokenizer.nextTokenAsDouble();
            currProp.cod[1] = tokenizer.nextTokenAsDouble();
            currProp.cod[2] = tokenizer.nextTokenAsDouble();

            currProp.Ddtt(0,0) = tokenizer.nextTokenAsDouble();
            currProp.Ddtt(0,1) = tokenizer.nextTokenAsDouble();
            currProp.Ddtt(0,2) = tokenizer.nextTokenAsDouble();
            currProp.Ddtt(1,0) = tokenizer.nextTokenAsDouble();
            currProp.Ddtt(1,1) = tokenizer.nextTokenAsDouble();
            currProp.Ddtt(1,2) = tokenizer.nextTokenAsDouble();
            currProp.Ddtt(2,0) = tokenizer.nextTokenAsDouble();
            currProp.Ddtt(2,1) = tokenizer.nextTokenAsDouble();
            currProp.Ddtt(2,2) = tokenizer.nextTokenAsDouble();

            currProp.Ddtr(0,0) = tokenizer.nextTokenAsDouble();
            currProp.Ddtr(0,1) = tokenizer.nextTokenAsDouble();
            currProp.Ddtr(0,2) = tokenizer.nextTokenAsDouble();
            currProp.Ddtr(1,0) = tokenizer.nextTokenAsDouble();
            currProp.Ddtr(1,1) = tokenizer.nextTokenAsDouble();
            currProp.Ddtr(1,2) = tokenizer.nextTokenAsDouble();
            currProp.Ddtr(2,0) = tokenizer.nextTokenAsDouble();
            currProp.Ddtr(2,1) = tokenizer.nextTokenAsDouble();
            currProp.Ddtr(2,2) = tokenizer.nextTokenAsDouble();

            currProp.Ddrr(0,0) = tokenizer.nextTokenAsDouble();
            currProp.Ddrr(0,1) = tokenizer.nextTokenAsDouble();
            currProp.Ddrr(0,2) = tokenizer.nextTokenAsDouble();
            currProp.Ddrr(1,0) = tokenizer.nextTokenAsDouble();
            currProp.Ddrr(1,1) = tokenizer.nextTokenAsDouble();
            currProp.Ddrr(1,2) = tokenizer.nextTokenAsDouble();
            currProp.Ddrr(2,0) = tokenizer.nextTokenAsDouble();
            currProp.Ddrr(2,1) = tokenizer.nextTokenAsDouble();
            currProp.Ddrr(2,2) = tokenizer.nextTokenAsDouble();                

            currProp.Xidtt(0,0) = tokenizer.nextTokenAsDouble();
            currProp.Xidtt(0,1) = tokenizer.nextTokenAsDouble();
            currProp.Xidtt(0,2) = tokenizer.nextTokenAsDouble();
            currProp.Xidtt(1,0) = tokenizer.nextTokenAsDouble();
            currProp.Xidtt(1,1) = tokenizer.nextTokenAsDouble();
            currProp.Xidtt(1,2) = tokenizer.nextTokenAsDouble();
            currProp.Xidtt(2,0) = tokenizer.nextTokenAsDouble();
            currProp.Xidtt(2,1) = tokenizer.nextTokenAsDouble();
            currProp.Xidtt(2,2) = tokenizer.nextTokenAsDouble();

            currProp.Xidrt(0,0) = tokenizer.nextTokenAsDouble();
            currProp.Xidrt(0,1) = tokenizer.nextTokenAsDouble();
            currProp.Xidrt(0,2) = tokenizer.nextTokenAsDouble();
            currProp.Xidrt(1,0) = tokenizer.nextTokenAsDouble();
            currProp.Xidrt(1,1) = tokenizer.nextTokenAsDouble();
            currProp.Xidrt(1,2) = tokenizer.nextTokenAsDouble();
            currProp.Xidrt(2,0) = tokenizer.nextTokenAsDouble();
            currProp.Xidrt(2,1) = tokenizer.nextTokenAsDouble();
            currProp.Xidrt(2,2) = tokenizer.nextTokenAsDouble();
            
            currProp.Xidtr(0,0) = tokenizer.nextTokenAsDouble();
            currProp.Xidtr(0,1) = tokenizer.nextTokenAsDouble();
            currProp.Xidtr(0,2) = tokenizer.nextTokenAsDouble();
            currProp.Xidtr(1,0) = tokenizer.nextTokenAsDouble();
            currProp.Xidtr(1,1) = tokenizer.nextTokenAsDouble();
            currProp.Xidtr(1,2) = tokenizer.nextTokenAsDouble();
            currProp.Xidtr(2,0) = tokenizer.nextTokenAsDouble();
            currProp.Xidtr(2,1) = tokenizer.nextTokenAsDouble();
            currProp.Xidtr(2,2) = tokenizer.nextTokenAsDouble();

            currProp.Xidrr(0,0) = tokenizer.nextTokenAsDouble();
            currProp.Xidrr(0,1) = tokenizer.nextTokenAsDouble();
            currProp.Xidrr(0,2) = tokenizer.nextTokenAsDouble();
            currProp.Xidrr(1,0) = tokenizer.nextTokenAsDouble();
            currProp.Xidrr(1,1) = tokenizer.nextTokenAsDouble();
            currProp.Xidrr(1,2) = tokenizer.nextTokenAsDouble();
            currProp.Xidrr(2,0) = tokenizer.nextTokenAsDouble();
            currProp.Xidrr(2,1) = tokenizer.nextTokenAsDouble();
            currProp.Xidrr(2,2) = tokenizer.nextTokenAsDouble(); 
            props.insert(std::map<std::string, HydroProp>::value_type(atomName, currProp));
        }
    }

    return props;
  }
  
  void LDForceManager::postCalculation() {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* integrableObject;
    Vector3d vel;
    Vector3d pos;
    Vector3d frc;
    Mat3x3d A;
    Mat3x3d Atrans;
    Vector3d Tb;
    Vector3d ji;
    double mass;
    unsigned int index = 0;
    for (mol = info_->beginMolecule(i); mol != NULL; mol = info_->nextMolecule(i)) {
      for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL;
	   integrableObject = mol->nextIntegrableObject(j)) {

          vel =integrableObject->getVel(); 
          if (integrableObject->isDirectional()){
             //calculate angular velocity in lab frame
             Mat3x3d I = integrableObject->getI();
             Vector3d angMom = integrableObject->getJ();
             Vector3d omega;

             if (integrableObject->isLinear()) {
                int linearAxis = integrableObject->linearAxis();
                int l = (linearAxis +1 )%3;
                int m = (linearAxis +2 )%3;
                omega[l] = angMom[l] /I(l, l);
                omega[m] = angMom[m] /I(m, m);
                
             } else {
                 omega[0] = angMom[0] /I(0, 0);
                 omega[1] = angMom[1] /I(1, 1);
                 omega[2] = angMom[2] /I(2, 2);
             }

             //apply friction force and torque at center of diffusion
             A = integrableObject->getA();
             Atrans = A.transpose();
             Vector3d rcd = Atrans * hydroProps_[index].cod;  
             Vector3d vcd = vel + cross(omega, rcd);
             vcd = A* vcd;
             Vector3d frictionForce = -(hydroProps_[index].Xidtt * vcd + hydroProps_[index].Xidrt * omega);
             frictionForce = Atrans*frictionForce;
             integrableObject->addFrc(frictionForce);
             Vector3d frictionTorque = - (hydroProps_[index].Xidtr * vcd + hydroProps_[index].Xidrr * omega);
             frictionTorque = Atrans*frictionTorque;
             integrableObject->addTrq(frictionTorque+ cross(rcd, frictionForce));
             
             //apply random force and torque at center of diffustion
             Vector3d randomForce;
             Vector3d randomTorque;
             genRandomForceAndTorque(randomForce, randomTorque, index, variance_);
             randomForce = Atrans*randomForce;
             randomTorque = Atrans* randomTorque;
             integrableObject->addFrc(randomForce);            
             integrableObject->addTrq(randomTorque + cross(rcd, randomForce ));
             
          } else {
             //spheric atom
             Vector3d frictionForce = -(hydroProps_[index].Xidtt *vel);     
             Vector3d randomForce;
             Vector3d randomTorque;
             genRandomForceAndTorque(randomForce, randomTorque, index, variance_);

             //randomForce /= OOPSEConstant::energyConvert;
             //randomTorque /= OOPSEConstant::energyConvert;
             integrableObject->addFrc(frictionForce+randomForce);             
          }

        ++index;
    
      }
    }    

    ForceManager::postCalculation();



  }

void LDForceManager::genRandomForceAndTorque(Vector3d& force, Vector3d& torque, unsigned int index, double variance) {
    /*
    SquareMatrix<double, 6> Dd;
    SquareMatrix<double, 6> S;
    Vector<double, 6> Z;
    Vector<double, 6> generalForce;
    Dd.setSubMatrix(0, 0, hydroProps_[index].Ddtt);
    Dd.setSubMatrix(0, 3, hydroProps_[index].Ddtr.transpose());
    Dd.setSubMatrix(3, 0, hydroProps_[index].Ddtr);
    Dd.setSubMatrix(3, 3, hydroProps_[index].Ddrr);
    CholeskyDecomposition(Dd, S);
    */

    SquareMatrix<double, 6> Xid;
    SquareMatrix<double, 6> S;
    Vector<double, 6> Z;
    Vector<double, 6> generalForce;
    Xid.setSubMatrix(0, 0, hydroProps_[index].Xidtt);
    Xid.setSubMatrix(0, 3, hydroProps_[index].Xidrt);
    Xid.setSubMatrix(3, 0, hydroProps_[index].Xidtr);
    Xid.setSubMatrix(3, 3, hydroProps_[index].Xidrr);
    CholeskyDecomposition(Xid, S);

    /*
    Xid *= variance;
    Z[0] = randNumGen_.randNorm(0, 1.0);
    Z[1] = randNumGen_.randNorm(0, 1.0);
    Z[2] = randNumGen_.randNorm(0, 1.0);
    Z[3] = randNumGen_.randNorm(0, 1.0);
    Z[4] = randNumGen_.randNorm(0, 1.0);
    Z[5] = randNumGen_.randNorm(0, 1.0);
    */
        
    Z[0] = randNumGen_.randNorm(0, variance);
    Z[1] = randNumGen_.randNorm(0, variance);
    Z[2] = randNumGen_.randNorm(0, variance);
    Z[3] = randNumGen_.randNorm(0, variance);
    Z[4] = randNumGen_.randNorm(0, variance);
    Z[5] = randNumGen_.randNorm(0, variance);
     

    generalForce = S*Z;
    
    force[0] = generalForce[0];
    force[1] = generalForce[1];
    force[2] = generalForce[2];
    torque[0] = generalForce[3];
    torque[1] = generalForce[4];
    torque[2] = generalForce[5];
    
}

}
