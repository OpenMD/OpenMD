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

#include "lattice/shapedLattice.hpp"

#include <cstdlib>

#include "brains/Register.hpp"
#include "lattice/LatticeFactory.hpp"
#include "math/Vector3.hpp"

using namespace std;
namespace OpenMD {
  shapedLattice::shapedLattice(RealType latticeConstant, string latticeType) :
      latticeConstant_(latticeConstant), latticeType_(latticeType) {
    registerLattice();
    simpleLattice_ = LatticeFactory::getInstance().createLattice(latticeType);
    if (simpleLattice_ == NULL) {
      std::cerr << "shapedLattice:: Error creating lattice" << std::endl;
      exit(1);
    }

    // Set the lattice constant
    std::vector<RealType> lc;
    lc.push_back(latticeConstant_);
    simpleLattice_->setLatticeConstant(lc);
    sitesComputed_ = false;
  }

  void shapedLattice::setGridDimension(Vector3d dimension) {
    dimension_ = dimension;
    // Find	number of unit cells in each direction
    beginNx_       = -(int)ceil(0.5 * dimension_[0] / latticeConstant_);
    beginNy_       = -(int)ceil(0.5 * dimension_[1] / latticeConstant_);
    beginNz_       = -(int)ceil(0.5 * dimension_[2] / latticeConstant_);
    endNx_         = (int)ceil(0.5 * dimension_[0] / latticeConstant_);
    endNy_         = (int)ceil(0.5 * dimension_[1] / latticeConstant_);
    endNz_         = (int)ceil(0.5 * dimension_[2] / latticeConstant_);
    sitesComputed_ = false;
  }

  void shapedLattice::setOrigin(Vector3d origin) {
    origin_ = origin;
    simpleLattice_->setOrigin(origin_);
    sitesComputed_ = false;
  }

  void shapedLattice::findSites() {
    sites_.clear();
    orientations_.clear();

    std::vector<Vector3d> latticePos;
    std::vector<Vector3d> pointsOrt = simpleLattice_->getLatticePointsOrt();
    int numMolPerCell               = simpleLattice_->getNumSitesPerCell();

    for (int i = beginNx_; i <= endNx_; i++) {
      for (int j = beginNy_; j <= endNy_; j++) {
        for (int k = beginNz_; k <= endNz_; k++) {
          // get the position of the cell sites
          simpleLattice_->getLatticePointsPos(latticePos, i, j, k);
          for (int l = 0; l < numMolPerCell; l++) {
            if (isInterior(latticePos[l])) {
              Vector3d myPoint = latticePos[l];
              Vector3d myOrt   = pointsOrt[l];
              sites_.push_back(myPoint);
              orientations_.push_back(myOrt);
            }
          }
        }
      }
    }
    sitesComputed_ = true;
  }

  std::vector<Vector3d> shapedLattice::getSites() {
    if (!sitesComputed_) { findSites(); }
    return sites_;
  }

  std::vector<Vector3d> shapedLattice::getOrientations() {
    if (!sitesComputed_) { findSites(); }
    return orientations_;
  }
}  // namespace OpenMD
