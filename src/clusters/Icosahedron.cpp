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

#include "clusters/Icosahedron.hpp"

#include <tuple>

namespace OpenMD {

  Icosahedron::Icosahedron() {
    Basis.clear();
    Edges.clear();
    Facets.clear();
    Points.clear();

    //
    // Initialize Basis vectors.
    //
    const RealType HGR = (sqrt(5.0) + 1.0) / 4.0;  // half of the golden ratio

    Basis.push_back(Vector3d(HGR, 0.0, 0.5));
    Basis.push_back(Vector3d(HGR, 0.0, -0.5));
    Basis.push_back(Vector3d(0.5, HGR, 0.0));
    Basis.push_back(Vector3d(-0.5, HGR, 0.0));
    Basis.push_back(Vector3d(0.0, 0.5, HGR));
    Basis.push_back(Vector3d(0.0, -0.5, HGR));
    Basis.push_back(Vector3d(0.5, -HGR, 0.0));
    Basis.push_back(Vector3d(0.0, 0.5, -HGR));
    Basis.push_back(Vector3d(-HGR, 0.0, 0.5));
    Basis.push_back(Vector3d(0.0, -0.5, -HGR));
    Basis.push_back(Vector3d(-HGR, 0.0, -0.5));
    Basis.push_back(Vector3d(-0.5, -HGR, 0.0));

    //
    // Initialize 30 edges
    //

    Edges.push_back(std::make_pair(0, 1));
    Edges.push_back(std::make_pair(0, 2));
    Edges.push_back(std::make_pair(0, 4));
    Edges.push_back(std::make_pair(0, 5));
    Edges.push_back(std::make_pair(0, 6));

    Edges.push_back(std::make_pair(10, 3));
    Edges.push_back(std::make_pair(10, 7));
    Edges.push_back(std::make_pair(10, 8));
    Edges.push_back(std::make_pair(10, 9));
    Edges.push_back(std::make_pair(10, 11));

    Edges.push_back(std::make_pair(1, 2));
    Edges.push_back(std::make_pair(1, 6));
    Edges.push_back(std::make_pair(1, 7));
    Edges.push_back(std::make_pair(1, 9));

    Edges.push_back(std::make_pair(8, 3));
    Edges.push_back(std::make_pair(8, 4));
    Edges.push_back(std::make_pair(8, 5));
    Edges.push_back(std::make_pair(8, 11));

    Edges.push_back(std::make_pair(2, 3));
    Edges.push_back(std::make_pair(2, 4));
    Edges.push_back(std::make_pair(2, 7));

    Edges.push_back(std::make_pair(11, 5));
    Edges.push_back(std::make_pair(11, 6));
    Edges.push_back(std::make_pair(11, 9));

    Edges.push_back(std::make_pair(6, 5));
    Edges.push_back(std::make_pair(6, 9));

    Edges.push_back(std::make_pair(3, 4));
    Edges.push_back(std::make_pair(3, 7));

    Edges.push_back(std::make_pair(7, 9));

    Edges.push_back(std::make_pair(5, 4));

    //
    // Initialize 20 facets
    //
    Facets.push_back(std::make_tuple(0, 1, 2));
    Facets.push_back(std::make_tuple(0, 2, 4));
    Facets.push_back(std::make_tuple(0, 4, 5));
    Facets.push_back(std::make_tuple(0, 5, 6));
    Facets.push_back(std::make_tuple(0, 1, 6));

    Facets.push_back(std::make_tuple(10, 3, 7));
    Facets.push_back(std::make_tuple(10, 3, 8));
    Facets.push_back(std::make_tuple(10, 8, 11));
    Facets.push_back(std::make_tuple(10, 9, 11));
    Facets.push_back(std::make_tuple(10, 7, 9));

    Facets.push_back(std::make_tuple(1, 2, 7));
    Facets.push_back(std::make_tuple(1, 7, 9));
    Facets.push_back(std::make_tuple(1, 6, 9));

    Facets.push_back(std::make_tuple(8, 5, 11));
    Facets.push_back(std::make_tuple(8, 4, 5));
    Facets.push_back(std::make_tuple(8, 3, 4));

    Facets.push_back(std::make_tuple(2, 3, 7));
    Facets.push_back(std::make_tuple(2, 3, 4));

    Facets.push_back(std::make_tuple(11, 5, 6));
    Facets.push_back(std::make_tuple(11, 6, 9));
  }

  int Icosahedron::getNpoints(int n) {
    int count = 0;
    for (int i = 0; i <= n; i++)
      count += np(i);
    return count;
  }

  int Icosahedron::np(int n) {
    if (n < 0)
      return -1;
    else if (n == 0)
      return 1;
    else if (n == 1)
      return 12;
    else if (n == 2)
      return 42;
    else {
      int count = 0;
      count += 12;            // edge particles
      count += (n - 1) * 30;  // side particles
      for (int i = 1; i <= n - 2; i++)
        count += i * 20;  // body particles
      return count;
    }
  }

  std::vector<Vector3d> Icosahedron::ih(int n) {
    if (n < 0) return Points;

    if (n == 0) {
      // center particle only

      Points.push_back(Vector3d(0.0, 0.0, 0.0));
      return Points;
    }

    //
    // Generate edge particles
    //
    for (std::vector<Vector3d>::iterator i = Basis.begin(); i != Basis.end();
         ++i) {
      Points.push_back((*i) * RealType(n));
    }

    //
    // Generate side particles
    //
    if (n < 2) return Points;

    for (std::vector<std::pair<int, int>>::iterator i = Edges.begin();
         i != Edges.end(); ++i) {
      Vector3d e1 = Basis[(*i).first] * RealType(n);
      Vector3d e2 = Basis[(*i).second] * RealType(n);

      for (int j = 1; j <= n - 1; j++) {
        Points.push_back(e1 + (e2 - e1) * RealType(j) / RealType(n));
      }
    }

    //
    // Generate body particles
    //
    if (n < 3) return Points;

    for (const auto& [first, second, third] : Facets) {
      Vector3d e1 = Basis[first] * RealType(n);
      Vector3d e2 = Basis[second] * RealType(n);
      Vector3d e3 = Basis[third] * RealType(n);

      for (int j = 1; j <= n - 2; j++) {
        Vector3d v1 = e1 + (e2 - e1) * RealType(j + 1) / RealType(n);
        Vector3d v2 = e1 + (e3 - e1) * RealType(j + 1) / RealType(n);

        for (int k = 1; k <= j; k++) {
          Points.push_back(v1 + (v2 - v1) * RealType(k) / RealType(j + 1));
        }
      }
    }

    return Points;
  }

  std::vector<Vector3d> Icosahedron::getPoints(int nshells) {
    // generate the coordinates
    for (int i = 0; i <= nshells; i++)
      ih(i);
    return Points;
  }
}  // namespace OpenMD
