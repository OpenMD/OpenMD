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
      /* Convert pressure from atm -> amu/(fs^2*Ang)*/
      targetPressure_ = simParams->getTargetPressure()/OOPSEConstant::pressureConvert;
    }

   
    if (simParams->getUsePeriodicBoundaryConditions()) {
      sprintf(painCave.errMsg, "SMIPDynamics error: You can't use the SMIPD integrator\n"
	      "   with periodic boundary conditions !\n");
      
      painCave.isFatal = 1;
      simError();
    } 

    if (!simParams->haveViscosity()) {
      sprintf(painCave.errMsg, "You can't use SMIPDynamics without a viscosity!\n");
      painCave.isFatal = 1;
      painCave.severity = OOPSE_ERROR;
      simError();
    }else{
      viscosity_ = simParams->getViscosity();
    }

    dt_ = simParams->getDt();

    

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


    Molecule* mol;
    StuntDouble* integrableObject;
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;		

    

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

    
    std::vector<Triangle> sMesh = surfaceMesh_->getMesh();
    int nTriangles = sMesh.size();



     /* Compute variance for random forces */
    std::vector<RealType>  randNums = genTriangleForces(nTriangles, 1.0);

    /* Loop over the mesh faces and apply random force to each of the faces*/    
    std::vector<Triangle>::iterator face;
    std::vector<StuntDouble*>::iterator vertex;
    int thisNumber = 0;
    for (face = sMesh.begin(); face != sMesh.end(); ++face){
     
      Triangle thisTriangle = *face;
      std::vector<StuntDouble*> vertexSDs = thisTriangle.getVertices();
      RealType thisArea = thisTriangle.getArea(); 
      // RealType sigma_t = sqrt(NumericConstant::PI/2.0)*((targetPressure_)*thisArea) /OOPSEConstant::energyConvert;
      // gamma_t_ = (NumericConstant::PI * targetPressure_* targetPressure_ * thisArea * thisArea * simParams->getDt()) /(4.0 * OOPSEConstant::kB*simParams->getTargetTemp());

      /* Get Random Force */
      Vector3d unitNormal = thisTriangle.getNormal();
      unitNormal.normalize();
      //Vector3d randomForce = -randNums[thisNumber] * sigma_t * unitNormal;
      Vector3d centroid = thisTriangle.getCentroid();
      Vector3d facetVel = thisTriangle.getFacetVelocity();
      RealType hydroLength = thisTriangle.getIncircleRadius()*2.0/NumericConstant::PI;

      RealType f_normal = viscosity_*hydroLength*OOPSEConstant::viscoConvert;
      RealType extPressure = -(targetPressure_ * thisArea)/OOPSEConstant::energyConvert;
      RealType randomForce = randNums[thisNumber++] * sqrt(2.0 * f_normal * OOPSEConstant::kb*targetTemp_/dt_);

      RealType dragForce = -f_normal * dot(facetVel, unitNormal);


      Vector3d langevinForce = (extPressure + randomForce + dragForce) * unitNormal;
      
      //      Vector3d dragForce = - gamma_t_ * dot(facetVel, unitNormal) * unitNormal / OOPSEConstant::energyConvert;
      
      // std::cout << " " << randomForce << " " << f_normal <<   std::endl;

      /* Apply triangle force to stuntdouble vertices */
      for (vertex = vertexSDs.begin(); vertex != vertexSDs.end(); ++vertex){
	if ((*vertex) != NULL){
	  // mass = integrableObject->getMass();
	  Vector3d vertexForce = langevinForce/3.0;
	  (*vertex)->addFrc(vertexForce);

	  if ((*vertex)->isDirectional()){

	    Vector3d vertexPos = (*vertex)->getPos();
	    Vector3d vertexCentroidVector = vertexPos - centroid;
	    (*vertex)->addTrq(cross(vertexCentroidVector,vertexForce));
	  }
	}  
      }
    } 

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
	//gaussRand[i] = fabs(randNumGen_.randNorm(0.0, 1.0));     
	gaussRand[i] = randNumGen_.randNorm(0.0, 1.0);     
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
