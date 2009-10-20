/*
 * Copyright (c) 2008, 2009 The University of Notre Dame. All Rights Reserved.
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
#include "utils/OOPSEConstant.hpp"
#include "math/ConvexHull.hpp"
#include "math/Triangle.hpp"

namespace oopse {

  SMIPDForceManager::SMIPDForceManager(SimInfo* info) : ForceManager(info) {

    simParams = info->getSimParams();
    veloMunge = new Velocitizer(info);
    
    // Create Hull, Convex Hull for now, other options later.
    
    surfaceMesh_ = new ConvexHull();
    
    /* Check that the simulation has target pressure and target
       temperature set */
    
    if (!simParams->haveTargetTemp()) {
      sprintf(painCave.errMsg, 
              "SMIPDynamics error: You can't use the SMIPD integrator\n"
	      "\twithout a targetTemp (K)!\n");      
      painCave.isFatal = 1;
      painCave.severity = OOPSE_ERROR;
      simError();
    } else {
      targetTemp_ = simParams->getTargetTemp();
    }
    
    if (!simParams->haveTargetPressure()) {
      sprintf(painCave.errMsg, 
              "SMIPDynamics error: You can't use the SMIPD integrator\n"
	      "\twithout a targetPressure (atm)!\n");      
      painCave.isFatal = 1;
      simError();
    } else {
      // Convert pressure from atm -> amu/(fs^2*Ang)
      targetPressure_ = simParams->getTargetPressure() / 
        OOPSEConstant::pressureConvert;
    }
   
    if (simParams->getUsePeriodicBoundaryConditions()) {
      sprintf(painCave.errMsg, 
              "SMIPDynamics error: You can't use the SMIPD integrator\n"
              "\twith periodic boundary conditions!\n");    
      painCave.isFatal = 1;
      simError();
    } 
    
    if (!simParams->haveThermalConductivity()) {
      sprintf(painCave.errMsg, 
              "SMIPDynamics error: You can't use the SMIPD integrator\n"
	      "\twithout a thermalConductivity (W m^-1 K^-1)!\n");
      painCave.isFatal = 1;
      painCave.severity = OOPSE_ERROR;
      simError();
    }else{
      thermalConductivity_ = simParams->getThermalConductivity() * 
        OOPSEConstant::thermalConductivityConvert;
    }

    if (!simParams->haveThermalLength()) {
      sprintf(painCave.errMsg, 
              "SMIPDynamics error: You can't use the SMIPD integrator\n"
	      "\twithout a thermalLength (Angstroms)!\n");
      painCave.isFatal = 1;
      painCave.severity = OOPSE_ERROR;
      simError();
    }else{
      thermalLength_ = simParams->getThermalLength();
    }
    
    dt_ = simParams->getDt();
  
    variance_ = 2.0 * OOPSEConstant::kb * targetTemp_ / dt_;

    // Build a vector of integrable objects to determine if the are
    // surface atoms
    Molecule* mol;
    StuntDouble* integrableObject;
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;

    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {          
      for (integrableObject = mol->beginIntegrableObject(j); 
           integrableObject != NULL;
           integrableObject = mol->nextIntegrableObject(j)) {	
	localSites_.push_back(integrableObject);
      }
    }   
  }  
   
  void SMIPDForceManager::postCalculation(bool needStress){
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* integrableObject;
  
    // Compute surface Mesh

    surfaceMesh_->computeHull(localSites_);

    // Get total area and number of surface stunt doubles
    RealType area = surfaceMesh_->getArea();
    std::vector<Triangle> sMesh = surfaceMesh_->getMesh();
    int nTriangles = sMesh.size();

    // Generate all of the necessary random forces
    std::vector<RealType>  randNums = genTriangleForces(nTriangles, variance_);

    // Loop over the mesh faces and apply external pressure to each 
    // of the faces
    std::vector<Triangle>::iterator face;
    std::vector<StuntDouble*>::iterator vertex;
    int thisFacet = 0;
    for (face = sMesh.begin(); face != sMesh.end(); ++face){
     
      Triangle thisTriangle = *face;
      std::vector<StuntDouble*> vertexSDs = thisTriangle.getVertices();
      RealType thisArea = thisTriangle.getArea(); 
      Vector3d unitNormal = thisTriangle.getNormal();
      unitNormal.normalize();
      Vector3d centroid = thisTriangle.getCentroid();
      Vector3d facetVel = thisTriangle.getFacetVelocity();
      RealType thisMass = thisTriangle.getFacetMass();

      // gamma is the drag coefficient normal to the face of the triangle      
      RealType gamma = thermalConductivity_ * thisMass * thisArea 
        / (2.0 * thermalLength_ * OOPSEConstant::kB);
      
      RealType extPressure = - (targetPressure_ * thisArea) / 
        OOPSEConstant::energyConvert;

      RealType randomForce = randNums[thisFacet++] * sqrt(gamma);
      RealType dragForce = -gamma * dot(facetVel, unitNormal);

      Vector3d langevinForce = (extPressure + randomForce + dragForce) * 
        unitNormal;
      
      // Apply triangle force to stuntdouble vertices
      for (vertex = vertexSDs.begin(); vertex != vertexSDs.end(); ++vertex){
	if ((*vertex) != NULL){
	  Vector3d vertexForce = langevinForce / 3.0;
          //          std::cout << "Adding force: " << facetVel << " to global id: " << (*vertex)->getGlobalIndex() << std::endl; 
	  (*vertex)->addFrc(vertexForce);	   
	}  
      }
    } 
    
    veloMunge->removeComDrift();
    veloMunge->removeAngularDrift();
    
    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    currSnapshot->setVolume(surfaceMesh_->getVolume());    
    ForceManager::postCalculation(needStress);   
  }
  
  
  std::vector<RealType> SMIPDForceManager::genTriangleForces(int nTriangles, 
                                                             RealType variance)
  {
    
    // zero fill the random vector before starting:
    std::vector<RealType> gaussRand;
    gaussRand.resize(nTriangles);
    std::fill(gaussRand.begin(), gaussRand.end(), 0.0);
    
#ifdef IS_MPI
    if (worldRank == 0) {
#endif
      for (int i = 0; i < nTriangles; i++) {
        gaussRand[i] = randNumGen_.randNorm(0.0, variance);
      }
#ifdef IS_MPI
    }
#endif
    
    // push these out to the other processors
    
#ifdef IS_MPI
    if (worldRank == 0) {
      MPI::COMM_WORLD.Bcast(&gaussRand[0], nTriangles, MPI::REALTYPE, 0);
    } else {
      MPI::COMM_WORLD.Bcast(&gaussRand[0], nTriangles, MPI::REALTYPE, 0);
    }
#endif
    
    return gaussRand;
  }
}
