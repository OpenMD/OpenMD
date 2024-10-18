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

#include "clusters/Cuboctahedron.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "utils/CaseConversion.hpp"

using namespace std;

namespace OpenMD {

  bool pairCompare(const pair<RealType, int>& l, const pair<RealType, int>& r) {
    return l.first < r.first;
  }

  Cuboctahedron::Cuboctahedron(std::string lattice, int cells, int planes) :
      lattice_(lattice), L_(cells), M_(planes) {
    Basis.clear();
    Points.clear();

    //
    // Initialize Basis vectors.
    //
    toUpper(lattice);

    if (lattice == "BCC") {
      Basis.push_back(Vector3d(0.0, 0.0, 0.0));
      Basis.push_back(Vector3d(0.5, 0.5, 0.5));
    } else if (lattice == "SC") {
      Basis.push_back(Vector3d(0.0, 0.0, 0.0));
    } else {
      // Default is FCC:
      Basis.push_back(Vector3d(0.0, 0.0, 0.0));
      Basis.push_back(Vector3d(0.5, 0.5, 0.0));
      Basis.push_back(Vector3d(0.5, 0.0, 0.5));
      Basis.push_back(Vector3d(0.0, 0.5, 0.5));
    }
  }

  vector<Vector3d> Cuboctahedron::getPoints() {
    // center of cluster

    Vector3d c(0.0);

    Vector3d d;
    vector<Vector3d> rawPoints;
    vector<pair<RealType, int>> dists;
    int idx;

    for (int i = -L_; i <= L_; i++) {
      for (int j = -L_; j <= L_; j++) {
        for (int k = -L_; k <= L_; k++) {
          for (vector<Vector3d>::iterator l = Basis.begin(); l != Basis.end();
               ++l) {
            Vector3d point = (*l) + Vector3d(i, j, k);
            if (inCluster(point)) {
              rawPoints.push_back(point);
              d   = point - c;
              idx = dists.size();
              dists.push_back(make_pair(d.lengthSquare(), idx));
            }
          }
        }
      }
    }

    // Sort the atoms using distance from center of cluster:

    sort(dists.begin(), dists.end(), pairCompare);

    for (vector<pair<RealType, int>>::iterator i = dists.begin();
         i != dists.end(); ++i) {
      Points.push_back(rawPoints[(*i).second] - c);
    }

    return Points;
  }

  bool Cuboctahedron::inCluster111(Vector3d r) {
    Vector3d c   = r.abs();
    RealType rad = 0.5 * RealType(L_);

    if ((c.x() < rad) && (c.y() < rad) && (c.z() < rad) &&
        (c.x() + c.y() + c.z() < RealType(M_)))
      return true;
    else
      return false;
  }

  bool Cuboctahedron::inCluster(Vector3d r) { return inCluster111(r); }
}  // namespace OpenMD
