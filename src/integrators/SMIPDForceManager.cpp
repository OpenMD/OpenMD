/*
 * Copyright (c) 2008 The University of Notre Dame. All Rights Reserved.
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
#include <iostream>
#include "integrators/SMIPDForceManager.hpp"
#include "math/CholeskyDecomposition.hpp"
#include "utils/OOPSEConstant.hpp"
#include "hydrodynamics/Sphere.hpp"
#include "hydrodynamics/Ellipsoid.hpp"
#include "utils/ElementsTable.hpp"
#include "math/ConvexHull.hpp"
#include "math/Triangle.hpp"


namespace oopse {

  SMIPDForceManager::SMIPDForceManager(SimInfo* info) : ForceManager(info), forceTolerance_(1e-6), maxIterNum_(4) {
    simParams = info->getSimParams();
    veloMunge = new Velocitizer(info);
    
    // Create Hull, Convex Hull for now, other options later.
    surfaceMesh_ = new ConvexHull();

    
    
    /* Check that the simulation has target pressure and target
       temperature set*/

    if (!simParams->haveTargetTemp()) {
      sprintf(painCave.errMsg, "You can't use the SMIPDynamics integrator without a targetTemp!\n");
      painCave.isFatal = 1;
      painCave.severity = OOPSE_ERROR;
      simError();
    } else {
      targetTemp_ = simParams->getTargetTemp();
    }

    if (!simParams->haveTargetPressure()) {
      sprintf(painCave.errMsg, "SMIPDynamics error: You can't use the SMIPD integrator\n"
	      "   without a targetPressure!\n");
      
      painCave.isFatal = 1;
      simError();
    } else {
      targetPressure_ = simParams->getTargetPressure()/OOPSEConstant::pressureConvert;
    }

   
    if (simParams->getUsePeriodicBoundaryConditions()) {
      sprintf(painCave.errMsg, "SMIPDynamics error: You can't use the SMIPD integrator\n"
	      "   with periodic boundary conditions !\n");
      
      painCave.isFatal = 1;
      simError();
    } 



    

    //Compute initial hull
    /*
    surfaceMesh_->computeHull(localSites_);
    Area0_ = surfaceMesh_->getArea();
    int nTriangles = surfaceMesh_->getNMeshElements();
    //    variance_ = 2.0 * OOPSEConstant::kb*simParams->getTargetTemp()/simParams->getDt();
    gamma_0_ = (NumericConstant::PI * targetPressure_* targetPressure_ * Area0_ * Area0_ * simParams->getDt()) /
      (4.0 * nTriangles * nTriangles* OOPSEConstant::kb*simParams->getTargetTemp());
    //RealType eta0 = gamma_0/
    */

    // Build the hydroProp map:
    std::map<std::string, HydroProp*> hydroPropMap;

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
                 "HydroPropFile must be set to a file name if SMIPDynamics\n"
                 "\tis specified for rigidBodies which contain more than one atom\n"
                 "\tTo create a HydroPropFile, run the \"Hydro\" program.\n\n"
		 "\tFor use with SMIPD, the default viscosity in Hydro should be\n"
		 "\tset to 1.0 because the friction and random forces will be\n"
                 "\tdynamically re-set assuming this is true.\n");
        painCave.severity = OOPSE_ERROR;
        painCave.isFatal = 1;
        simError();  
      }      

      for (mol = info->beginMolecule(i); mol != NULL; 
           mol = info->nextMolecule(i)) {
        for (integrableObject = mol->beginIntegrableObject(j); 
             integrableObject != NULL;
             integrableObject = mol->nextIntegrableObject(j)) {

          std::map<std::string, HydroProp*>::iterator iter = hydroPropMap.find(integrableObject->getType());
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
      
      std::map<std::string, HydroProp*> hydroPropMap;
      for (mol = info->beginMolecule(i); mol != NULL; 
           mol = info->nextMolecule(i)) {
        for (integrableObject = mol->beginIntegrableObject(j); 
             integrableObject != NULL;
             integrableObject = mol->nextIntegrableObject(j)) {
          Shape* currShape = NULL;

          if (integrableObject->isAtom()){
            Atom* atom = static_cast<Atom*>(integrableObject);
            AtomType* atomType = atom->getAtomType();
            if (atomType->isGayBerne()) {
              DirectionalAtomType* dAtomType = dynamic_cast<DirectionalAtomType*>(atomType);              
              GenericData* data = dAtomType->getPropertyByName("GayBerne");
              if (data != NULL) {
                GayBerneParamGenericData* gayBerneData = dynamic_cast<GayBerneParamGenericData*>(data);
                
                if (gayBerneData != NULL) {  
                  GayBerneParam gayBerneParam = gayBerneData->getData();
                  currShape = new Ellipsoid(V3Zero, 
                                            gayBerneParam.GB_l / 2.0, 
                                            gayBerneParam.GB_d / 2.0, 
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
            } else {
              if (atomType->isLennardJones()){
                GenericData* data = atomType->getPropertyByName("LennardJones");
                if (data != NULL) {
                  LJParamGenericData* ljData = dynamic_cast<LJParamGenericData*>(data);
                  if (ljData != NULL) {
                    LJParam ljParam = ljData->getData();
                    currShape = new Sphere(atom->getPos(), 2.0);
                  } else {
                    sprintf( painCave.errMsg,
                             "Can not cast GenericData to LJParam\n");
                    painCave.severity = OOPSE_ERROR;
                    painCave.isFatal = 1;
                    simError();          
                  }       
                }
              } else {
                int aNum = etab.GetAtomicNum((atom->getType()).c_str());
                if (aNum != 0) {
                  currShape = new Sphere(atom->getPos(), 2.0);
                } else {
                  sprintf( painCave.errMsg,
                           "Could not find atom type in default element.txt\n");
                  painCave.severity = OOPSE_ERROR;
                  painCave.isFatal = 1;
                  simError();          
                }
              }
            }
          }
          HydroProp* currHydroProp = currShape->getHydroProp(1.0,simParams->getTargetTemp());
          std::map<std::string, HydroProp*>::iterator iter = hydroPropMap.find(integrableObject->getType());
          if (iter != hydroPropMap.end()) 
            hydroProps_.push_back(iter->second);
          else {
            currHydroProp->complete();
            hydroPropMap.insert(std::map<std::string, HydroProp*>::value_type(integrableObject->getType(), currHydroProp));
            hydroProps_.push_back(currHydroProp);
          }
        }
      }
    }

    /* Compute hull first time through to get the area of t=0*/

    //Build a vector of integrable objects to determine if the are surface atoms 
    for (mol = info_->beginMolecule(i); mol != NULL; mol = info_->nextMolecule(i)) {          
      for (integrableObject = mol->beginIntegrableObject(j); integrableObject != NULL;
           integrableObject = mol->nextIntegrableObject(j)) {	
	localSites_.push_back(integrableObject);
      }
    }


  }  

  std::map<std::string, HydroProp*> SMIPDForceManager::parseFrictionFile(const std::string& filename) {
    std::map<std::string, HydroProp*> props;
    std::ifstream ifs(filename.c_str());
    if (ifs.is_open()) {
      
    }
    
    const unsigned int BufferSize = 65535;
    char buffer[BufferSize];   
    while (ifs.getline(buffer, BufferSize)) {
      HydroProp* currProp = new HydroProp(buffer);
      props.insert(std::map<std::string, HydroProp*>::value_type(currProp->getName(), currProp));
    }

    return props;
  }
   
  void SMIPDForceManager::postCalculation(bool needStress){
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* integrableObject;
    RealType mass;
    Vector3d pos;
    Vector3d frc;
    Mat3x3d A;
    Mat3x3d Atrans;
    Vector3d Tb;
    Vector3d ji;
    unsigned int index = 0;
    int fdf;
   
    fdf = 0;
  
    /*Compute surface Mesh*/
    surfaceMesh_->computeHull(localSites_);

    /* Get area and number of surface stunt doubles and compute new variance */
     RealType area = surfaceMesh_->getArea();
     int nSurfaceSDs = surfaceMesh_->getNs();

    
    std::vector<Triangle*> sMesh = surfaceMesh_->getMesh();
    int nTriangles = sMesh.size();



     /* Compute variance for random forces */
   
    RealType sigma_t = sqrt(NumericConstant::PI/2.0)*((targetPressure_)*area/nTriangles)
       /OOPSEConstant::energyConvert;

    gamma_t_ = (NumericConstant::PI * targetPressure_* targetPressure_ * area * area * simParams->getDt()) /
      (4.0 * nTriangles * nTriangles* OOPSEConstant::kB*simParams->getTargetTemp()) /OOPSEConstant::energyConvert;

    std::vector<RealType>  randNums = genTriangleForces(nTriangles, sigma_t);

    /* Loop over the mesh faces and apply random force to each of the faces*/
    
    
    std::vector<Triangle*>::iterator face;
    std::vector<StuntDouble*>::iterator vertex;
    int thisNumber = 0;
    for (face = sMesh.begin(); face != sMesh.end(); ++face){
     
      Triangle* thisTriangle = *face;
      std::vector<StuntDouble*> vertexSDs = thisTriangle->getVertices();
      
      /* Get Random Force */
      Vector3d unitNormal = thisTriangle->getNormal();
      unitNormal.normalize();
      Vector3d randomForce = -randNums[thisNumber] * unitNormal;
      Vector3d centroid = thisTriangle->getCentroid();

      Vector3d langevinForce = randomForce - gamma_t_*thisTriangle->getFacetVelocity();
      
      for (vertex = vertexSDs.begin(); vertex != vertexSDs.end(); ++vertex){
	if ((*vertex) != NULL){
	  // mass = integrableObject->getMass();
	  Vector3d vertexForce = langevinForce/3.0;
	  (*vertex)->addFrc(vertexForce);
	  if (integrableObject->isDirectional()){
	    Vector3d vertexPos = (*vertex)->getPos();
	    Vector3d vertexCentroidVector = vertexPos - centroid;
	    (*vertex)->addTrq(cross(vertexCentroidVector,vertexForce));
	  }
	}  
      }
    } 

    /* Now loop over all surface particles and apply the drag*/
    /*
    std::vector<StuntDouble*> surfaceSDs = surfaceMesh_->getSurfaceAtoms();
    for (vertex = surfaceSDs.begin(); vertex != surfaceSDs.end(); ++vertex){
      integrableObject = *vertex;
      mass = integrableObject->getMass();
      if (integrableObject->isDirectional()){
	
	// preliminaries for directional objects:
	
	A = integrableObject->getA();
	Atrans = A.transpose();
	Vector3d rcrLab = Atrans * hydroProps_[index]->getCOR();  
	//apply random force and torque at center of resistance
	Mat3x3d I = integrableObject->getI();
	Vector3d omegaBody;
	
	// What remains contains velocity explicitly, but the velocity required
	// is at the full step: v(t + h), while we have initially the velocity
	// at the half step: v(t + h/2).  We need to iterate to converge the
	// friction force and friction torque vectors.
	
	// this is the velocity at the half-step:
        
	Vector3d vel =integrableObject->getVel();
	Vector3d angMom = integrableObject->getJ();
	
	//estimate velocity at full-step using everything but friction forces:           
	
	frc = integrableObject->getFrc();
	Vector3d velStep = vel + (dt2_ /mass * OOPSEConstant::energyConvert) * frc;
	
	Tb = integrableObject->lab2Body(integrableObject->getTrq());
	Vector3d angMomStep = angMom + (dt2_ * OOPSEConstant::energyConvert) * Tb;                             
	
	Vector3d omegaLab;
	Vector3d vcdLab;
	Vector3d vcdBody;
	Vector3d frictionForceBody;
	Vector3d frictionForceLab(0.0);
	Vector3d oldFFL;  // used to test for convergence
	Vector3d frictionTorqueBody(0.0);
	Vector3d oldFTB;  // used to test for convergence
	Vector3d frictionTorqueLab;
	RealType fdot;
	RealType tdot;

	//iteration starts here:
	
	for (int k = 0; k < maxIterNum_; k++) {
	  
	  if (integrableObject->isLinear()) {
	    int linearAxis = integrableObject->linearAxis();
	    int l = (linearAxis +1 )%3;
	    int m = (linearAxis +2 )%3;
	    omegaBody[l] = angMomStep[l] /I(l, l);
	    omegaBody[m] = angMomStep[m] /I(m, m);
            
	  } else {
	    omegaBody[0] = angMomStep[0] /I(0, 0);
	    omegaBody[1] = angMomStep[1] /I(1, 1);
	    omegaBody[2] = angMomStep[2] /I(2, 2);
	  }
          
	  omegaLab = Atrans * omegaBody;
          
	  // apply friction force and torque at center of resistance
          
	  vcdLab = velStep + cross(omegaLab, rcrLab);       
	  vcdBody = A * vcdLab;
	  frictionForceBody = -(hydroProps_[index]->getXitt() * vcdBody + hydroProps_[index]->getXirt() * omegaBody);
	  oldFFL = frictionForceLab;
	  frictionForceLab = Atrans * frictionForceBody;
	  oldFTB = frictionTorqueBody;
	  frictionTorqueBody = -(hydroProps_[index]->getXitr() * vcdBody + hydroProps_[index]->getXirr() * omegaBody);
	  frictionTorqueLab = Atrans * frictionTorqueBody;
          
	  // re-estimate velocities at full-step using friction forces:
              
	  velStep = vel + (dt2_ / mass * OOPSEConstant::energyConvert) * (frc + frictionForceLab);
	  angMomStep = angMom + (dt2_ * OOPSEConstant::energyConvert) * (Tb + frictionTorqueBody);

	  // check for convergence (if the vectors have converged, fdot and tdot will both be 1.0):
              
	  fdot = dot(frictionForceLab, oldFFL) / frictionForceLab.lengthSquare();
	  tdot = dot(frictionTorqueBody, oldFTB) / frictionTorqueBody.lengthSquare();
	  
	  if (fabs(1.0 - fdot) <= forceTolerance_ && fabs(1.0 - tdot) <= forceTolerance_)
	    break; // iteration ends here
	}
	
	integrableObject->addFrc(frictionForceLab);
	integrableObject->addTrq(frictionTorqueLab + cross(rcrLab, frictionForceLab));

            
      } else {
 	//spherical atom

	// What remains contains velocity explicitly, but the velocity required
	// is at the full step: v(t + h), while we have initially the velocity
	// at the half step: v(t + h/2).  We need to iterate to converge the
	// friction force vector.
	
	// this is the velocity at the half-step:
        
	Vector3d vel =integrableObject->getVel();
	
	//estimate velocity at full-step using everything but friction forces:           
	
	frc = integrableObject->getFrc();
	Vector3d velStep = vel + (dt2_ / mass * OOPSEConstant::energyConvert) * frc;
	
	Vector3d frictionForce(0.0);
	Vector3d oldFF;  // used to test for convergence
	RealType fdot;
	
	//iteration starts here:
	
	for (int k = 0; k < maxIterNum_; k++) {
	  
	  oldFF = frictionForce;                            
	  frictionForce = -hydroProps_[index]->getXitt() * velStep;
	  //frictionForce = -gamma_t*velStep;
	  // re-estimate velocities at full-step using friction forces:
          
	  velStep = vel + (dt2_ / mass * OOPSEConstant::energyConvert) * (frc + frictionForce);
	  
	  // check for convergence (if the vector has converged, fdot will be 1.0):
          
	  fdot = dot(frictionForce, oldFF) / frictionForce.lengthSquare();
          
	  if (fabs(1.0 - fdot) <= forceTolerance_)
	    break; // iteration ends here
	}
	
	integrableObject->addFrc(frictionForce);
	
        
      }
  
      
  }
    */
    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    currSnapshot->setVolume(surfaceMesh_->getVolume());
    
    ForceManager::postCalculation(needStress);   
  }

  void SMIPDForceManager::genRandomForceAndTorque(Vector3d& force, Vector3d& torque, unsigned int index, RealType variance) {

    
    Vector<RealType, 6> Z;
    Vector<RealType, 6> generalForce;
    

    Z[0] = randNumGen_.randNorm(0, variance);
    Z[1] = randNumGen_.randNorm(0, variance);
    Z[2] = randNumGen_.randNorm(0, variance);
    Z[3] = randNumGen_.randNorm(0, variance);
    Z[4] = randNumGen_.randNorm(0, variance);
    Z[5] = randNumGen_.randNorm(0, variance);
     
    generalForce = hydroProps_[index]->getS()*Z;
    
    force[0] = generalForce[0];
    force[1] = generalForce[1];
    force[2] = generalForce[2];
    torque[0] = generalForce[3];
    torque[1] = generalForce[4];
    torque[2] = generalForce[5];
    
} 
  std::vector<RealType> SMIPDForceManager::genTriangleForces(int nTriangles, RealType variance) {

    // zero fill the random vector before starting:
    std::vector<RealType> gaussRand;
    gaussRand.resize(nTriangles);
    std::fill(gaussRand.begin(), gaussRand.end(), 0.0);


#ifdef IS_MPI
    if (worldRank == 0) {
#endif
      for (int i = 0; i < nTriangles; i++) {
	gaussRand[i] = fabs(randNumGen_.randNorm(0.0, 1.0));     
      }
#ifdef IS_MPI
    }
#endif

    // push these out to the other processors

#ifdef IS_MPI
    if (worldRank == 0) {
      MPI_Bcast(&gaussRand[0], nTriangles, MPI_REAL, 0, MPI_COMM_WORLD);
    } else {
      MPI_Bcast(&gaussRand[0], nTriangles, MPI_REAL, 0, MPI_COMM_WORLD);
    }
#endif

    for (int i = 0; i < nTriangles; i++) {
      gaussRand[i] = gaussRand[i] * variance;
    }

    return gaussRand;
  }





}
