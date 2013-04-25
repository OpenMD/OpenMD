/*
 * Copyright (c) 2008, 2009, 2010 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
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
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *                                                                      
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
#include <fstream> 
#include <iostream>
#include "integrators/LangevinHullForceManager.hpp"
#include "utils/PhysicalConstants.hpp"
#include "math/ConvexHull.hpp"
#include "math/AlphaHull.hpp"
#include "math/Triangle.hpp"
#include "math/CholeskyDecomposition.hpp"
#ifdef IS_MPI
#include <mpi.h>
#endif

using namespace std;
namespace OpenMD {
  
  LangevinHullForceManager::LangevinHullForceManager(SimInfo* info) : 
    ForceManager(info) {
   
    simParams = info->getSimParams();
    veloMunge = new Velocitizer(info);
    
    // Create Hull, Convex Hull for now, other options later.
    
    stringToEnumMap_["Convex"] = hullConvex;
    stringToEnumMap_["AlphaShape"] = hullAlphaShape;
    stringToEnumMap_["Unknown"] = hullUnknown;
    
    const std::string ht = simParams->getHULL_Method();
    
    std::map<std::string, HullTypeEnum>::iterator iter;
    iter = stringToEnumMap_.find(ht);
    hullType_ = (iter == stringToEnumMap_.end()) ? 
      LangevinHullForceManager::hullUnknown : iter->second;

    switch(hullType_) {
    case hullConvex :
      surfaceMesh_ = new ConvexHull();
      break;
    case hullAlphaShape :
      surfaceMesh_ = new AlphaHull(simParams->getAlpha());
      break;
    case hullUnknown :
    default :
      sprintf(painCave.errMsg, 
              "LangevinHallForceManager: Unknown Hull_Method was requested!\n");
      painCave.isFatal = 1;
      simError();      
      break;
    }

    doThermalCoupling_ = true;
    doPressureCoupling_ = true;

    /* Check that the simulation has target pressure and target
       temperature set */    
    if (!simParams->haveTargetTemp()) {
      sprintf(painCave.errMsg, 
              "LangevinHullForceManager: no targetTemp (K) was set.\n"
              "\tOpenMD is turning off the thermal coupling to the bath.\n");
      painCave.isFatal = 0;
      painCave.severity = OPENMD_INFO;
      simError();
      doThermalCoupling_ = false;
    } else {
      targetTemp_ = simParams->getTargetTemp();

      if (!simParams->haveViscosity()) {
        sprintf(painCave.errMsg, 
                "LangevinHullForceManager: no viscosity was set.\n"
                "\tOpenMD is turning off the thermal coupling to the bath.\n");
        painCave.isFatal = 0;
        painCave.severity = OPENMD_INFO;
        simError();      
        doThermalCoupling_ = false;
      }else{
        viscosity_ = simParams->getViscosity();
        if ( fabs(viscosity_) < 1e-6 ) {
          sprintf(painCave.errMsg, 
                  "LangevinHullDynamics: The bath viscosity was set lower\n"
                  "\tthan 1e-6 poise.  OpenMD is turning off the thermal\n"
                  "\tcoupling to the bath.\n");
          painCave.isFatal = 0;
          painCave.severity = OPENMD_INFO;
          simError();
          doThermalCoupling_ = false;
        }      
      }      
    }
    if (!simParams->haveTargetPressure()) {
      sprintf(painCave.errMsg, 
              "LangevinHullForceManager: no targetPressure (atm) was set.\n"
              "\tOpenMD is turning off the pressure coupling to the bath.\n");
      painCave.isFatal = 0;
      painCave.severity = OPENMD_INFO;
      simError();
      doPressureCoupling_ = false;
    } else {
      // Convert pressure from atm -> amu/(fs^2*Ang)
      targetPressure_ = simParams->getTargetPressure() / 
        PhysicalConstants::pressureConvert;
    }
    if (simParams->getUsePeriodicBoundaryConditions()) {
      sprintf(painCave.errMsg, 
              "LangevinHallForceManager: You can't use the Langevin Hull\n"
              "\tintegrator with periodic boundary conditions turned on!\n");
      painCave.isFatal = 1;
      simError();
    } 

    dt_ = simParams->getDt();

    if (doThermalCoupling_)
      variance_ = 2.0 * PhysicalConstants::kb * targetTemp_ / dt_;

    // Build a vector of integrable objects to determine if the are
    // surface atoms
    Molecule* mol;
    StuntDouble* sd;
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    
    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {          
      for (sd = mol->beginIntegrableObject(j); 
           sd != NULL;
           sd = mol->nextIntegrableObject(j)) {	
	localSites_.push_back(sd);
      }
    }
    
    // We need to make an initial guess at the bounding box in order
    // to compute long range forces in ForceMatrixDecomposition:

    // Compute surface Mesh
    surfaceMesh_->computeHull(localSites_);
    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
  }  
  
  void LangevinHullForceManager::postCalculation(){
  
    int nTriangles, thisFacet;
    RealType area, thisArea, thisMass;
    vector<Triangle> sMesh;
    Triangle thisTriangle;
    vector<Triangle>::iterator face;
    vector<StuntDouble*> vertexSDs;
    vector<StuntDouble*>::iterator vertex;

    Vector3d unitNormal, centroid, facetVel;
    Vector3d langevinForce, vertexForce;
    Vector3d extPressure, randomForce, dragForce;

    Mat3x3d hydroTensor, S;
    vector<Vector3d> randNums;

    // Compute surface Mesh
    surfaceMesh_->computeHull(localSites_);
    // Get total area and number of surface stunt doubles
    area = surfaceMesh_->getArea();
    sMesh = surfaceMesh_->getMesh();
    nTriangles = sMesh.size();

    if (doThermalCoupling_) {
      // Generate all of the necessary random forces
      randNums = genTriangleForces(nTriangles, variance_);
    }
    
    // Loop over the mesh faces
    thisFacet = 0;
    for (face = sMesh.begin(); face != sMesh.end(); ++face){
      thisTriangle = *face;
      vertexSDs = thisTriangle.getVertices();
      thisArea = thisTriangle.getArea(); 
      unitNormal = thisTriangle.getUnitNormal();
      centroid = thisTriangle.getCentroid();
      facetVel = thisTriangle.getFacetVelocity();
      thisMass = thisTriangle.getFacetMass();

      langevinForce = V3Zero;

      if (doPressureCoupling_) {
        extPressure = -unitNormal * (targetPressure_ * thisArea) /
          PhysicalConstants::energyConvert;
        langevinForce += extPressure;
      }

      if (doThermalCoupling_) {
        hydroTensor = thisTriangle.computeHydrodynamicTensor(viscosity_);      
        hydroTensor *= PhysicalConstants::viscoConvert;
        CholeskyDecomposition(hydroTensor, S);
        randomForce = S * randNums[thisFacet++];
        dragForce = -hydroTensor * facetVel;
        langevinForce += randomForce + dragForce;
      }
      
      // Apply triangle force to stuntdouble vertices
      for (vertex = vertexSDs.begin(); vertex != vertexSDs.end(); ++vertex){
	if ((*vertex) != NULL){
          vertexForce = langevinForce / RealType(3.0);
	  (*vertex)->addFrc(vertexForce);	   
	}  
      }
    } 
    
    veloMunge->removeComDrift();
    veloMunge->removeAngularDrift();
    
    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    currSnapshot->setVolume(surfaceMesh_->getVolume());   
    currSnapshot->setHullVolume(surfaceMesh_->getVolume());
    ForceManager::postCalculation();   
  }
    
  vector<Vector3d> LangevinHullForceManager::genTriangleForces(int nTriangles, 
                                                               RealType var) {
    // zero fill the random vector before starting:
    vector<Vector3d> gaussRand;
    gaussRand.resize(nTriangles);
    std::fill(gaussRand.begin(), gaussRand.end(), V3Zero);
    
#ifdef IS_MPI
    if (worldRank == 0) {
#endif
      for (int i = 0; i < nTriangles; i++) {
        gaussRand[i][0] = randNumGen_.randNorm(0.0, var);
        gaussRand[i][1] = randNumGen_.randNorm(0.0, var);
        gaussRand[i][2] = randNumGen_.randNorm(0.0, var);
      }
#ifdef IS_MPI
    }
#endif
    
    // push these out to the other processors
    
#ifdef IS_MPI
    if (worldRank == 0) {
      MPI::COMM_WORLD.Bcast(&gaussRand[0], nTriangles*3, MPI::REALTYPE, 0);
    } else {
      MPI::COMM_WORLD.Bcast(&gaussRand[0], nTriangles*3, MPI::REALTYPE, 0);
    }
#endif
    
    return gaussRand;
  }
}
