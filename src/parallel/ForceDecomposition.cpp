/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "parallel/ForceDecomposition.hpp"
#include "math/Vector3.hpp"

using namespace std;
namespace OpenMD {

  ForceDecomposition::ForceDecomposition(SimInfo* info,
                                         InteractionManager* iMan) :
    info_(info), interactionMan_(iMan), needVelocities_(false) {

    sman_ = info_->getSnapshotManager();
    storageLayout_ = sman_->getStorageLayout();
    ff_ = info_->getForceField();

    usePeriodicBoundaryConditions_ = info->getSimParams()->getUsePeriodicBoundaryConditions();

    Globals* simParams_ = info_->getSimParams();    
    if (simParams_->havePrintHeatFlux()) {
      if (simParams_->getPrintHeatFlux()) {
        needVelocities_ = true;
      }
    }

    if (simParams_->haveSkinThickness()) {
      skinThickness_ = simParams_->getSkinThickness();
    } else {      
      skinThickness_ = 1.0;
      sprintf(painCave.errMsg,
              "ForceDecomposition: No value was set for the skinThickness.\n"
              "\tOpenMD will use a default value of %f Angstroms\n"
              "\tfor this simulation\n", skinThickness_);
      painCave.severity = OPENMD_INFO;
      painCave.isFatal = 0;
      simError();
    }             

    // cellOffsets are the partial space for the cell lists used in
    // constructing the neighbor lists
    cellOffsets_.clear();
    cellOffsets_.push_back( Vector3i(0, 0, 0) );
    cellOffsets_.push_back( Vector3i(1, 0, 0) );
    cellOffsets_.push_back( Vector3i(1, 1, 0) );
    cellOffsets_.push_back( Vector3i(0, 1, 0) );
    cellOffsets_.push_back( Vector3i(-1,1, 0) );
    cellOffsets_.push_back( Vector3i(0, 0, 1) );
    cellOffsets_.push_back( Vector3i(1, 0, 1) );
    cellOffsets_.push_back( Vector3i(1, 1, 1) );
    cellOffsets_.push_back( Vector3i(0, 1, 1) );
    cellOffsets_.push_back( Vector3i(-1,1, 1) );
    cellOffsets_.push_back( Vector3i(-1,0, 1) );
    cellOffsets_.push_back( Vector3i(-1,-1,1) );
    cellOffsets_.push_back( Vector3i(0, -1,1) );
    cellOffsets_.push_back( Vector3i(1, -1,1) );
  }

  void ForceDecomposition::setCutoffRadius(RealType rcut) {
    rCut_ = rcut;
    rList_ = rCut_ + skinThickness_;
    rListSq_ = rList_ * rList_;
  }

  void ForceDecomposition::fillPreForceData(SelfData &sdat, int atom) {

    sdat.atid = idents[atom];
    sdat.selfPot = selfPot;
    sdat.selePot = selectedSelfPot;
    
    if (storageLayout_ & DataStorage::dslDensity) {
      sdat.rho = snap_->atomData.density[atom];
    }
    
    if (storageLayout_ & DataStorage::dslParticlePot) {
      sdat.particlePot = snap_->atomData.particlePot[atom];
    }
    
  }

  void ForceDecomposition::fillSelfData(SelfData &sdat, int atom) {

    sdat.atid = idents[atom];
    sdat.selfPot = selfPot;
    sdat.selePot = selectedSelfPot;
    
    if (storageLayout_ & DataStorage::dslSkippedCharge) {
      sdat.skippedCharge = snap_->atomData.skippedCharge[atom];
    }

    if (storageLayout_ & DataStorage::dslParticlePot) {
      sdat.particlePot = snap_->atomData.particlePot[atom];
    }
    
    if (storageLayout_ & DataStorage::dslFlucQPosition) {
      sdat.flucQ = snap_->atomData.flucQPos[atom];
    }
    
    if (storageLayout_ & DataStorage::dslFlucQForce) {
      sdat.flucQfrc = snap_->atomData.flucQFrc[atom];
    }
  }
  
  void ForceDecomposition::unpackPreForceData(SelfData &sdat, int atom) {

    selfPot = sdat.selfPot;
    selectedSelfPot = sdat.selePot;
    
    if (storageLayout_ & DataStorage::dslFunctional) {
      snap_->atomData.functional[atom] += sdat.frho;
    }
    
    if (storageLayout_ & DataStorage::dslFunctionalDerivative) {
      snap_->atomData.functionalDerivative[atom] += sdat.dfrhodrho;
    }

    if (sdat.doParticlePot && (storageLayout_ & DataStorage::dslParticlePot)) {
      snap_->atomData.particlePot[atom] = sdat.particlePot;
    }
    
  }

  void ForceDecomposition::unpackSelfData(SelfData &sdat, int atom) {

    selfPot = sdat.selfPot;
    selectedSelfPot = sdat.selePot;
    
    if (storageLayout_ & DataStorage::dslFlucQForce) {
      snap_->atomData.flucQFrc[atom] = sdat.flucQfrc;
    }    
  }

  bool ForceDecomposition::checkNeighborList() {
    RealType st2 = pow( skinThickness_ / 2.0, 2);
    std::size_t nGroups = snap_->cgData.position.size();
    if (needVelocities_) 
      snap_->cgData.setStorageLayout(DataStorage::dslPosition |
                                     DataStorage::dslVelocity);
    
    // if we have changed the group identities or haven't set up the
    // saved positions we automatically will need a neighbor list update:
    
    if ( saved_CG_positions_.size() != nGroups ) return true;

    RealType dispmax = 0.0;
    Vector3d disp;    

    for (std::size_t i = 0; i < nGroups; i++) {
      disp = snap_->cgData.position[i]  - saved_CG_positions_[i];
      dispmax = max(dispmax, disp.lengthSquare());
    }

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &dispmax, 1, MPI_REALTYPE, MPI_MAX, 
                  MPI_COMM_WORLD);
#endif

    return (dispmax > st2) ? true : false;
  }

  void ForceDecomposition::addToHeatFlux(Vector3d hf) {
    Vector3d chf = snap_->getConductiveHeatFlux();
    chf += hf;
    snap_->setConductiveHeatFlux(chf);
  }
  void ForceDecomposition::setHeatFlux(Vector3d hf) {
    snap_->setConductiveHeatFlux(hf);
  }
}
