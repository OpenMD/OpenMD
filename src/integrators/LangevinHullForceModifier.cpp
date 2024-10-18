/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#include "integrators/LangevinHullForceModifier.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <string>

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "brains/ForceModifier.hpp"
#include "math/AlphaHull.hpp"
#include "math/CholeskyDecomposition.hpp"
#include "math/ConvexHull.hpp"
#include "math/Triangle.hpp"
#include "utils/Constants.hpp"

using namespace std;

namespace OpenMD {

  LangevinHullForceModifier::LangevinHullForceModifier(SimInfo* info) :
      ForceModifier {info} {
    simParams_ = info->getSimParams();
    veloMunge  = std::make_unique<Velocitizer>(info);

    // Create Hull, Convex Hull for now, other options later.

    stringToEnumMap_["Convex"]     = hullConvex;
    stringToEnumMap_["AlphaShape"] = hullAlphaShape;
    stringToEnumMap_["Unknown"]    = hullUnknown;

    const std::string ht = simParams_->getHULL_Method();

    std::map<std::string, HullTypeEnum>::iterator iter;
    iter      = stringToEnumMap_.find(ht);
    hullType_ = (iter == stringToEnumMap_.end()) ?
                    LangevinHullForceModifier::hullUnknown :
                    iter->second;

    switch (hullType_) {
    case hullConvex:
      surfaceMesh_ = new ConvexHull();
      break;
    case hullAlphaShape:
      surfaceMesh_ = new AlphaHull(simParams_->getAlpha());
      break;
    case hullUnknown:
    default:
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "LangevinHullForceModifier: Unknown Hull_Method was requested!\n");
      painCave.isFatal = 1;
      simError();
      break;
    }

    doThermalCoupling_  = true;
    doPressureCoupling_ = true;

    /* Check that the simulation has target pressure and target
     temperature set */
    if (!simParams_->haveTargetTemp()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "LangevinHullForceModifier: no targetTemp (K) was set.\n"
               "\tOpenMD is turning off the thermal coupling to the bath.\n");
      painCave.isFatal  = 0;
      painCave.severity = OPENMD_INFO;
      simError();
      doThermalCoupling_ = false;
    } else {
      targetTemp_ = simParams_->getTargetTemp();

      if (!simParams_->haveViscosity()) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "LangevinHullForceModifier: no viscosity was set.\n"
                 "\tOpenMD is turning off the thermal coupling to the bath.\n");
        painCave.isFatal  = 0;
        painCave.severity = OPENMD_INFO;
        simError();
        doThermalCoupling_ = false;
      } else {
        viscosity_ = simParams_->getViscosity();
        if (fabs(viscosity_) < 1e-6) {
          snprintf(
              painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
              "LangevinHullForceModifier: The bath viscosity was set lower\n"
              "\tthan 1e-6 poise.  OpenMD is turning off the thermal\n"
              "\tcoupling to the bath.\n");
          painCave.isFatal  = 0;
          painCave.severity = OPENMD_INFO;
          simError();
          doThermalCoupling_ = false;
        }
      }
    }
    if (!simParams_->haveTargetPressure()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "LangevinHullForceModifier: no targetPressure (atm) was set.\n"
               "\tOpenMD is turning off the pressure coupling to the bath.\n");
      painCave.isFatal  = 0;
      painCave.severity = OPENMD_INFO;
      simError();
      doPressureCoupling_ = false;
    } else {
      // Convert pressure from atm -> amu/(fs^2*Ang)
      targetPressure_ =
          simParams_->getTargetPressure() / Constants::pressureConvert;
    }
    if (simParams_->getUsePeriodicBoundaryConditions()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "LangevinHullForceModifier: You can't use the Langevin Hull\n"
               "\tintegrator with periodic boundary conditions turned on!\n");
      painCave.isFatal = 1;
      simError();
    }

    dt_ = simParams_->getDt();

#ifdef IS_MPI
    if (worldRank == 0) {
#endif
      if (doThermalCoupling_) {
        randNumGen_ = info->getRandomNumberGenerator();

        RealType stdDev = std::sqrt(2.0 * Constants::kb * targetTemp_ / dt_);

        forceDistribution_ = std::normal_distribution<RealType>(0.0, stdDev);
      }
#ifdef IS_MPI
    }
#endif

    // Build a vector of integrable objects to determine if the are
    // surface atoms
    Molecule* mol;
    StuntDouble* sd;
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        localSites_.push_back(sd);
      }
    }

    // We need to make an initial guess at the bounding box in order
    // to compute long range forces in ForceMatrixDecomposition:

    // Compute surface Mesh
    surfaceMesh_->computeHull(localSites_);
  }

  LangevinHullForceModifier::~LangevinHullForceModifier() {
    delete surfaceMesh_;
  }

  void LangevinHullForceModifier::modifyForces() {
    int nTriangles, thisFacet;
    RealType thisArea;
    std::vector<Triangle> sMesh;
    Triangle thisTriangle;
    std::vector<Triangle>::iterator face;
    std::vector<StuntDouble*> vertexSDs;
    std::vector<StuntDouble*>::iterator vertex;

    Vector3d unitNormal, centroid, facetVel;
    Vector3d langevinForce, vertexForce;
    Vector3d extPressure, randomForce, dragForce;

    Mat3x3d hydroTensor, S;
    std::vector<Vector3d> randNums;

    // Compute surface Mesh
    surfaceMesh_->computeHull(localSites_);
    // Get number of surface stunt doubles
    sMesh      = surfaceMesh_->getMesh();
    nTriangles = sMesh.size();

    if (doThermalCoupling_) {
      // Generate all of the necessary random forces
      randNums = genTriangleForces(nTriangles);
    }

    // Loop over the mesh faces
    thisFacet = 0;
    for (face = sMesh.begin(); face != sMesh.end(); ++face) {
      thisTriangle = *face;
      vertexSDs    = thisTriangle.getVertices();
      thisArea     = thisTriangle.getArea();
      unitNormal   = thisTriangle.getUnitNormal();
      centroid     = thisTriangle.getCentroid();
      facetVel     = thisTriangle.getFacetVelocity();

      langevinForce = V3Zero;

      if (doPressureCoupling_) {
        extPressure = -unitNormal * (targetPressure_ * thisArea) /
                      Constants::energyConvert;
        langevinForce += extPressure;
      }

      if (doThermalCoupling_) {
        hydroTensor = thisTriangle.computeHydrodynamicTensor(viscosity_);
        hydroTensor *= Constants::viscoConvert;
        CholeskyDecomposition(hydroTensor, S);
        randomForce = S * randNums[thisFacet++];
        dragForce   = -hydroTensor * facetVel;
        langevinForce += randomForce + dragForce;
      }

      // Apply triangle force to stuntdouble vertices
      for (vertex = vertexSDs.begin(); vertex != vertexSDs.end(); ++vertex) {
        if ((*vertex) != NULL) {
          vertexForce = langevinForce / RealType(3.0);
          (*vertex)->addFrc(vertexForce);
        }
      }
    }

    if (simParams_->getConserveLinearMomentum()) veloMunge->removeComDrift();
    if (simParams_->getConserveAngularMomentum())
      veloMunge->removeAngularDrift();

    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    currSnapshot->setVolume(surfaceMesh_->getVolume());
    currSnapshot->setHullVolume(surfaceMesh_->getVolume());
  }

  std::vector<Vector3d> LangevinHullForceModifier::genTriangleForces(
      int nTriangles) {
    // zero fill the random vector before starting:
    std::vector<Vector3d> gaussRand(nTriangles);
    std::fill(gaussRand.begin(), gaussRand.end(), V3Zero);

#ifdef IS_MPI
    if (worldRank == 0) {
#endif
      for (int i = 0; i < nTriangles; i++) {
        gaussRand[i][0] = forceDistribution_(*randNumGen_);
        gaussRand[i][1] = forceDistribution_(*randNumGen_);
        gaussRand[i][2] = forceDistribution_(*randNumGen_);
      }
#ifdef IS_MPI
    }
    // push these out to the other processors
    // Same command on all nodes:
    MPI_Bcast(gaussRand[0].getArrayPointer(), nTriangles * 3, MPI_REALTYPE, 0,
              MPI_COMM_WORLD);
#endif

    return gaussRand;
  }
}  // namespace OpenMD
