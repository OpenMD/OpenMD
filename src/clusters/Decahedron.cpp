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

#include "clusters/Decahedron.hpp"

#include <algorithm>
#include <cmath>

#include "utils/Constants.hpp"

using namespace std;

namespace OpenMD {

  Decahedron::Decahedron(int columnAtoms, int shells, int twinAtoms) :
      N_(columnAtoms), M_(shells), K_(twinAtoms) {
    Basis.clear();
    Points.clear();

    //
    // Initialize Basis vectors.
    //
    const RealType phi  = 2.0 * Constants::PI / 5.0;  // 72 degrees
    const RealType r3o2 = 0.5 * sqrt(3.0);

    Basis.push_back(
        Vector3d(r3o2 * cos(0.0 * phi), r3o2 * sin(0.0 * phi), 0.0));
    Basis.push_back(
        Vector3d(r3o2 * cos(1.0 * phi), r3o2 * sin(1.0 * phi), 0.0));
    Basis.push_back(
        Vector3d(r3o2 * cos(2.0 * phi), r3o2 * sin(2.0 * phi), 0.0));
    Basis.push_back(
        Vector3d(r3o2 * cos(3.0 * phi), r3o2 * sin(3.0 * phi), 0.0));
    Basis.push_back(
        Vector3d(r3o2 * cos(4.0 * phi), r3o2 * sin(4.0 * phi), 0.0));
  }

  vector<Vector3d> Decahedron::getPoints() {
    // Generate central column of Decahedron

    for (int i = 0; i < N_; i++) {
      Points.push_back(Vector3d(0.0, 0.0, RealType(i) - 0.5 * (N_ - 1)));
    }

    for (int i = 1; i < M_ + 1; i++) {
      // generate the shells of the decahedron:

      vector<Vector3d> ring;

      if (i > K_ - 1) {
        ring = truncatedRing(i, i - K_ + 1);
      } else {
        ring = truncatedRing(i, 0);
      }

      // shift the rings in the z-direction (along the shell)

      for (int j = 0; j < N_ - i; j++) {
        Vector3d shift =
            Vector3d(0, 0, -0.5 * RealType((N_ - i) - 1) + RealType(j));

        for (vector<Vector3d>::iterator k = ring.begin(); k != ring.end();
             ++k) {
          Points.push_back((*k) + shift);
        }
      }
    }
    return Points;
  }

  vector<Vector3d> Decahedron::truncatedRing(int n, int k) {
    // This function generates the rings of a Decahedron
    // n: index of shell (order of ring)
    // k: how many atoms are missing from both ends of one side of
    //    pentagon ring

    vector<Vector3d> ring;

    // Generate atomic coordinates along each side of pentagonal ring
    for (int i = 0; i < 5; i++) {
      Vector3d b1 = Basis[i];
      Vector3d b2 = Basis[(i + 1) % 5];

      if (k == 0) {
        // without truncation
        for (int j = 0; j < n; j++) {
          ring.push_back(RealType(n) * b1 + RealType(j) * (b2 - b1));
        }

      } else {
        for (int j = k; j <= n - k; j++) {
          // with truncation
          ring.push_back(RealType(n) * b1 + RealType(j) * (b2 - b1));
        }
      }
    }
    return ring;
  }

  CurlingStoneDecahedron::CurlingStoneDecahedron(int columnAtoms, int shells,
                                                 int twinAtoms,
                                                 int truncatedPlanes) :
      Decahedron(columnAtoms, shells, twinAtoms),
      T_(truncatedPlanes) {}

  vector<Vector3d> CurlingStoneDecahedron::getPoints() {
    vector<Vector3d> raw = Decahedron::getPoints();
    vector<Vector3d> snipped;
    RealType maxZ, minZ;

    maxZ = raw.begin()->z();
    minZ = raw.begin()->z();

    for (vector<Vector3d>::iterator i = raw.begin(); i != raw.end(); ++i) {
      maxZ = max(maxZ, (*i).z());
      minZ = min(minZ, (*i).z());
    }

    for (vector<Vector3d>::iterator i = raw.begin(); i != raw.end(); ++i) {
      if (((*i).z() < maxZ - 0.995 * (T_ / 2.0)) &&
          ((*i).z() > minZ + 0.995 * (T_ / 2.0))) {
        snipped.push_back((*i));
      }
    }
    return snipped;
  }

}  // namespace OpenMD
