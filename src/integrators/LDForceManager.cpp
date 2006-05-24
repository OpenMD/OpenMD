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
#include "hydrodynamics/Sphere.hpp"
#include "hydrodynamics/Ellipsoid.hpp"
#include "openbabel/mol.hpp"

using namespace OpenBabel;
namespace oopse {

  LDForceManager::LDForceManager(SimInfo* info) : ForceManager(info){
    Globals* simParams = info->getSimParams();
        
    sphericalBoundaryConditions_ = false;
    if (simParams->getUseSphericalBoundaryConditions()) {
      sphericalBoundaryConditions_ = true;
      if (simParams->haveLangevinBufferRadius()) {
        langevinBufferRadius_ = simParams->getLangevinBufferRadius();
      } else {
        sprintf( painCave.errMsg,
                 "langevinBufferRadius must be specified " 
                 "when useSphericalBoundaryConditions is turned on.\n");
        painCave.severity = OOPSE_ERROR;
        painCave.isFatal = 1;
        simError();  
      }
    
      if (simParams->haveFrozenBufferRadius()) {
        frozenBufferRadius_ = simParams->getFrozenBufferRadius();
      } else {
        sprintf( painCave.errMsg,
                 "frozenBufferRadius must be specified " 
                 "when useSphericalBoundaryConditions is turned on.\n");
        painCave.severity = OOPSE_ERROR;
        painCave.isFatal = 1;
        simError();  
      }

      if (frozenBufferRadius_ < langevinBufferRadius_) {
        sprintf( painCave.errMsg,
                 "frozenBufferRadius has been set smaller than the " 
                 "langevinBufferRadius.  This is probably an error.\n");
        painCave.severity = OOPSE_WARNING;
        painCave.isFatal = 0;
        simError();  
      }
    }

    // Build the hydroProp map:
    std::map<std::string, HydroProp> hydroPropMap;

    Molecule* mol;
    StuntDouble* integrableObject;
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;		
    bool needHydroPropFile = false;
    
    for (mol = info->beginMolecule(i); mol != NULL; 
         mol = info->nextMolecule(i)) {
      for (integrableObject = mol->beginIntegrableObject(j); 
           integrableObject != NULL;
           integrableObject = mol->nextIntegrableObject(j)) {
        
        if (integrableObject->isRigidBody()) {
          RigidBody* rb = static_cast<RigidBody*>(integrableObject);
          if (rb->getNumAtoms() > 1) needHydroPropFile = true;
        }
        
      }
    }
        

    if (needHydroPropFile) {               
      if (simParams->haveHydroPropFile()) {
        hydroPropMap = parseFrictionFile(simParams->getHydroPropFile());
      } else {               
        sprintf( painCave.errMsg,
                 "HydroPropFile must be set to a file name if Langevin\n"
                 "\tDynamics is specified for rigidBodies which contain more\n"
                 "\tthan one atom.  To create a HydroPropFile, run \"Hydro\".\n");
        painCave.severity = OOPSE_ERROR;
        painCave.isFatal = 1;
        simError();  
      }      

      for (mol = info->beginMolecule(i); mol != NULL; 
           mol = info->nextMolecule(i)) {
        for (integrableObject = mol->beginIntegrableObject(j); 
             integrableObject != NULL;
             integrableObject = mol->nextIntegrableObject(j)) {

          std::map<std::string, HydroProp>::iterator iter = hydroPropMap.find(integrableObject->getType());
          if (iter != hydroPropMap.end()) {
            hydroProps_.push_back(iter->second);
          } else {
            sprintf( painCave.errMsg,
                     "Can not find resistance tensor for atom [%s]\n", integrableObject->getType().c_str());
            painCave.severity = OOPSE_ERROR;
            painCave.isFatal = 1;
            simError();  
          }        
        }
      }
    } else {

      std::map<std::string, HydroProp> hydroPropMap;
      for (mol = info->beginMolecule(i); mol != NULL; 
           mol = info->nextMolecule(i)) {
        for (integrableObject = mol->beginIntegrableObject(j); 
             integrableObject != NULL;
             integrableObject = mol->nextIntegrableObject(j)) {
          Shape* currShape = NULL;
          if (integrableObject->isDirectionalAtom()) {
            DirectionalAtom* dAtom = static_cast<DirectionalAtom*>(integrableObject);
            AtomType* atomType = dAtom->getAtomType();
            if (atomType->isGayBerne()) {
              DirectionalAtomType* dAtomType = dynamic_cast<DirectionalAtomType*>(atomType);
              
              GenericData* data = dAtomType->getPropertyByName("GayBerne");
              if (data != NULL) {
                GayBerneParamGenericData* gayBerneData = dynamic_cast<GayBerneParamGenericData*>(data);
                
                if (gayBerneData != NULL) {  
                  GayBerneParam gayBerneParam = gayBerneData->getData();
                  currShape = new Ellipsoid(V3Zero, 
                                            gayBerneParam.GB_sigma/2.0, 
                                            gayBerneParam.GB_l2b_ratio*gayBerneParam.GB_sigma/2.0, 
                                            Mat3x3d::identity());
                } else {
                  sprintf( painCave.errMsg,
                           "Can not cast GenericData to GayBerneParam\n");
                  painCave.severity = OOPSE_ERROR;
                  painCave.isFatal = 1;
                  simError();   
                }
              } else {
                sprintf( painCave.errMsg, "Can not find Parameters for GayBerne\n");
                painCave.severity = OOPSE_ERROR;
                painCave.isFatal = 1;
                simError();    
              }
            }
          } else {
            Atom* atom = static_cast<Atom*>(integrableObject);
            AtomType* atomType = atom->getAtomType();
            if (atomType->isLennardJones()){
              GenericData* data = atomType->getPropertyByName("LennardJones");
              if (data != NULL) {
                LJParamGenericData* ljData = dynamic_cast<LJParamGenericData*>(data);
                
                if (ljData != NULL) {
                  LJParam ljParam = ljData->getData();
                  currShape = new Sphere(atom->getPos(), ljParam.sigma/2.0);
                } else {
                  sprintf( painCave.errMsg,
                           "Can not cast GenericData to LJParam\n");
                  painCave.severity = OOPSE_ERROR;
                  painCave.isFatal = 1;
                  simError();          
                }       
              }
            } else {
              int obanum = etab.GetAtomicNum((atom->getType()).c_str());
              if (obanum != 0) {
                currShape = new Sphere(atom->getPos(), etab.GetVdwRad(obanum));
              } else {
                sprintf( painCave.errMsg,
                         "Could not find atom type in default element.txt\n");
                painCave.severity = OOPSE_ERROR;
                painCave.isFatal = 1;
                simError();          
              }
            }
          }
          HydroProps currHydroProp = currShape->getHydroProps(simParams->getViscosity(),simParams->getTargetTemp());
          std::map<std::string, HydroProp>::iterator iter = hydroPropMap.find(integrableObject->getType());
          if (iter != hydroPropMap.end()) 
            hydroProps_.push_back(iter->second);
          else {
            HydroProp myProp;
            myProp.cor = V3Zero;
            for (int i1 = 0; i1 < 3; i1++) {
              for (int j1 = 0; j1 < 3; j1++) {
                myProp.Xirtt(i1,j1) = currHydroProp.Xi(i1,j1);
                myProp.Xirrt(i1,j1) = currHydroProp.Xi(i1,j1+3);
                myProp.Xirtr(i1,j1) = currHydroProp.Xi(i1+3,j1);
                myProp.Xirrr(i1,j1) = currHydroProp.Xi(i1+3,j1+3);
              }
            }
            CholeskyDecomposition(currHydroProp.Xi, myProp.S);
            hydroPropMap.insert(std::map<std::string, HydroProp>::value_type(integrableObject->getType(), myProp));
            hydroProps_.push_back(myProp);
          }
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
      if (tokenizer.countTokens() >= 40) {
        std::string atomName = tokenizer.nextToken();
        currProp.cor[0] = tokenizer.nextTokenAsDouble();
        currProp.cor[1] = tokenizer.nextTokenAsDouble();
        currProp.cor[2] = tokenizer.nextTokenAsDouble();
        
        currProp.Xirtt(0,0) = tokenizer.nextTokenAsDouble();
        currProp.Xirtt(0,1) = tokenizer.nextTokenAsDouble();
        currProp.Xirtt(0,2) = tokenizer.nextTokenAsDouble();
        currProp.Xirtt(1,0) = tokenizer.nextTokenAsDouble();
        currProp.Xirtt(1,1) = tokenizer.nextTokenAsDouble();
        currProp.Xirtt(1,2) = tokenizer.nextTokenAsDouble();
        currProp.Xirtt(2,0) = tokenizer.nextTokenAsDouble();
        currProp.Xirtt(2,1) = tokenizer.nextTokenAsDouble();
        currProp.Xirtt(2,2) = tokenizer.nextTokenAsDouble();
        
        currProp.Xirrt(0,0) = tokenizer.nextTokenAsDouble();
        currProp.Xirrt(0,1) = tokenizer.nextTokenAsDouble();
        currProp.Xirrt(0,2) = tokenizer.nextTokenAsDouble();
        currProp.Xirrt(1,0) = tokenizer.nextTokenAsDouble();
        currProp.Xirrt(1,1) = tokenizer.nextTokenAsDouble();
        currProp.Xirrt(1,2) = tokenizer.nextTokenAsDouble();
        currProp.Xirrt(2,0) = tokenizer.nextTokenAsDouble();
        currProp.Xirrt(2,1) = tokenizer.nextTokenAsDouble();
        currProp.Xirrt(2,2) = tokenizer.nextTokenAsDouble();
        
        currProp.Xirtr(0,0) = tokenizer.nextTokenAsDouble();
        currProp.Xirtr(0,1) = tokenizer.nextTokenAsDouble();
        currProp.Xirtr(0,2) = tokenizer.nextTokenAsDouble();
        currProp.Xirtr(1,0) = tokenizer.nextTokenAsDouble();
        currProp.Xirtr(1,1) = tokenizer.nextTokenAsDouble();
        currProp.Xirtr(1,2) = tokenizer.nextTokenAsDouble();
        currProp.Xirtr(2,0) = tokenizer.nextTokenAsDouble();
        currProp.Xirtr(2,1) = tokenizer.nextTokenAsDouble();
        currProp.Xirtr(2,2) = tokenizer.nextTokenAsDouble();
        
        currProp.Xirrr(0,0) = tokenizer.nextTokenAsDouble();
        currProp.Xirrr(0,1) = tokenizer.nextTokenAsDouble();
        currProp.Xirrr(0,2) = tokenizer.nextTokenAsDouble();
        currProp.Xirrr(1,0) = tokenizer.nextTokenAsDouble();
        currProp.Xirrr(1,1) = tokenizer.nextTokenAsDouble();
        currProp.Xirrr(1,2) = tokenizer.nextTokenAsDouble();
        currProp.Xirrr(2,0) = tokenizer.nextTokenAsDouble();
        currProp.Xirrr(2,1) = tokenizer.nextTokenAsDouble();
        currProp.Xirrr(2,2) = tokenizer.nextTokenAsDouble(); 
        
        SquareMatrix<RealType, 6> Xir;
        Xir.setSubMatrix(0, 0, currProp.Xirtt);
        Xir.setSubMatrix(0, 3, currProp.Xirrt);
        Xir.setSubMatrix(3, 0, currProp.Xirtr);
        Xir.setSubMatrix(3, 3, currProp.Xirrr);
        CholeskyDecomposition(Xir, currProp.S);            
        
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
    RealType mass;
    unsigned int index = 0;
    bool doLangevinForces;
    bool freezeMolecule;
    int fdf;
    
    fdf = 0;
    for (mol = info_->beginMolecule(i); mol != NULL; mol = info_->nextMolecule(i)) {

      doLangevinForces = true;           
      freezeMolecule = false;

      if (sphericalBoundaryConditions_) {
        
        Vector3d molPos = mol->getCom();
        RealType molRad = molPos.length();
        
        doLangevinForces = false;
        
        if (molRad > langevinBufferRadius_) { 
          doLangevinForces = true;
          freezeMolecule = false;
        }
        if (molRad > frozenBufferRadius_) {
          doLangevinForces = false;
          freezeMolecule = true;
        }
      }
      
      for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL;
           integrableObject = mol->nextIntegrableObject(j)) {
          
        if (freezeMolecule) 
          fdf += integrableObject->freeze();
        
        if (doLangevinForces) {          
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
            
            //apply friction force and torque at center of resistance
            A = integrableObject->getA();
            Atrans = A.transpose();
            Vector3d rcr = Atrans * hydroProps_[index].cor;  
            Vector3d vcdLab = vel + cross(omega, rcr);
            Vector3d vcdBody = A* vcdLab;
            Vector3d frictionForceBody = -(hydroProps_[index].Xirtt * vcdBody + hydroProps_[index].Xirrt * omega);
            Vector3d frictionForceLab = Atrans*frictionForceBody;
            integrableObject->addFrc(frictionForceLab);
            Vector3d frictionTorqueBody = - (hydroProps_[index].Xirtr * vcdBody + hydroProps_[index].Xirrr * omega);
            Vector3d frictionTorqueLab = Atrans*frictionTorqueBody;
            integrableObject->addTrq(frictionTorqueLab+ cross(rcr, frictionForceLab));
            
            //apply random force and torque at center of resistance
            Vector3d randomForceBody;
            Vector3d randomTorqueBody;
            genRandomForceAndTorque(randomForceBody, randomTorqueBody, index, variance_);
            Vector3d randomForceLab = Atrans*randomForceBody;
            Vector3d randomTorqueLab = Atrans* randomTorqueBody;
            integrableObject->addFrc(randomForceLab);            
            integrableObject->addTrq(randomTorqueLab + cross(rcr, randomForceLab ));             
            
          } else {
            //spherical atom
            Vector3d frictionForce = -(hydroProps_[index].Xirtt *vel);     
            Vector3d randomForce;
            Vector3d randomTorque;
            genRandomForceAndTorque(randomForce, randomTorque, index, variance_);
            
            integrableObject->addFrc(frictionForce+randomForce);             
          }
        }
          
        ++index;
    
      }
    }    
    info_->setFdf(fdf);
    
    ForceManager::postCalculation();   
  }

void LDForceManager::genRandomForceAndTorque(Vector3d& force, Vector3d& torque, unsigned int index, RealType variance) {


    Vector<RealType, 6> Z;
    Vector<RealType, 6> generalForce;

        
    Z[0] = randNumGen_.randNorm(0, variance);
    Z[1] = randNumGen_.randNorm(0, variance);
    Z[2] = randNumGen_.randNorm(0, variance);
    Z[3] = randNumGen_.randNorm(0, variance);
    Z[4] = randNumGen_.randNorm(0, variance);
    Z[5] = randNumGen_.randNorm(0, variance);
     

    generalForce = hydroProps_[index].S*Z;
    
    force[0] = generalForce[0];
    force[1] = generalForce[1];
    force[2] = generalForce[2];
    torque[0] = generalForce[3];
    torque[1] = generalForce[4];
    torque[2] = generalForce[5];
    
}

}
