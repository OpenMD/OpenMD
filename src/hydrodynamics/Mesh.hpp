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

#ifndef HYDRODYNAMICS_MESH_HPP_
#define HYDRODYNAMICS_MESH_HPP_

#include <config.h>

#include <cassert>
#include <string>
#include <vector>

#include "hydrodynamics/Shape.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Triangle.hpp"
#include "math/Vector3.hpp"

namespace OpenMD {
  class Mesh : public Shape {
  public:
    Mesh() {};
    virtual ~Mesh() {};
    // punt on this one for now:
    virtual bool isInterior(Vector3d pos) { return false; }
    virtual std::pair<Vector3d, Vector3d> getBoundingBox() {
      Vector3d bMax = facets_[0].vertex1();
      Vector3d bMin = facets_[0].vertex1();

      for (auto& t : facets_) {
        for (int i = 0; i < 3; i++) {
          bMax[i] = max(bMax[i], t.vertex1()[i]);
          bMin[i] = min(bMin[i], t.vertex1()[i]);
          bMax[i] = max(bMax[i], t.vertex2()[i]);
          bMin[i] = min(bMin[i], t.vertex2()[i]);
          bMax[i] = max(bMax[i], t.vertex3()[i]);
          bMin[i] = min(bMin[i], t.vertex3()[i]);
        }
      }
      return make_pair(bMax, bMin);
    }
    virtual bool hasAnalyticalSolution() { return false; }
    virtual HydroProp* getHydroProp(RealType viscosity) {
      HydroProp* props = new HydroProp();
      props->setCenterOfResistance(V3Zero);
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Mesh was asked to return an analytic HydroProps.\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
      return props;
    }
    virtual Vector3d getOrigin() { return Vector3d(0.0); }
    virtual bool isComposite() { return false; }
    virtual bool isSpherical() { return false; }
    virtual bool isMesh() { return true; }

    void add(Triangle t) { facets_.push_back(t); }

    void add(const Vector3d& vertex1, const Vector3d& vertex2,
             const Vector3d& vertex3) {
      Triangle t(vertex1, vertex2, vertex3);
      facets_.push_back(t);
    }

    void addQuadrilateral(const Vector3d& p1, const Vector3d& p2,
                          const Vector3d& p3, const Vector3d& p4) {
      add(p1, p2, p4);
      add(p4, p2, p3);
    };

    void clear() { facets_.clear(); }

    size_t size() { return facets_.size(); }

    Mesh& scale(const Vector3d& s) {
      double x, y, z;
      x = ((s.x() != 0) ? s.x() : 1);
      y = ((s.y() != 0) ? s.y() : 1);
      z = ((s.z() != 0) ? s.z() : 1);
      Vector3d _s(x, y, z);

      Vector3d v1, v2, v3;
      for (auto& t : facets_) {
        v1.Vmul(t.vertex1(), _s);
        v2.Vmul(t.vertex2(), _s);
        v3.Vmul(t.vertex3(), _s);
        t = {v1, v2, v3};
      }
      return *this;
    }

    Mesh& translate(const Vector3d& trans) {
      Vector3d v1, v2, v3;
      for (auto& t : facets_) {
        v1 = t.vertex1() + trans;
        v2 = t.vertex2() + trans;
        v3 = t.vertex3() + trans;
        t  = {v1, v2, v3};
      }
      return *this;
    }

    Mesh& flipNormal() {
      for (auto& t : facets_) {
        t.flipNormal();
      }
      return *this;
    }

    void stlWrite(const std::string& filename, const std::string& name = "") {
      std::ofstream file(filename);

      file << "solid"
           << " " << name << std::endl;
      for (auto& t : facets_) {
        file << "\t"
             << "facet normal"
             << " " << t.getUnitNormal() << std::endl;
        file << "\t\t"
             << "outer loop" << std::endl;
        file << "\t\t\t"
             << " "
             << "vertex"
             << " " << t.vertex1() << std::endl;
        file << "\t\t\t"
             << " "
             << "vertex"
             << " " << t.vertex2() << std::endl;
        file << "\t\t\t"
             << " "
             << "vertex"
             << " " << t.vertex3() << std::endl;
        file << "\t\t"
             << "endloop" << std::endl;
        file << "\t"
             << "endfacet" << std::endl;
      }
      file << "endsolid"
           << " " << name << std::endl;

      file.close();
    }

    Mesh& operator+=(const Mesh& m) {
      facets_.insert(facets_.end(), (m.facets_).begin(), (m.facets_).end());
      return *this;
    }

    inline void add(const Mesh& m1, const Mesh& m2) {
      this->facets_.insert(this->facets_.end(), (m1.facets_).begin(),
                           (m1.facets_).end());
      this->facets_.insert(this->facets_.end(), (m2.facets_).begin(),
                           (m2.facets_).end());
    }

    std::vector<Triangle> getFacets() { return facets_; }

  private:
    std::vector<Triangle> facets_;
  };

  inline Mesh operator+(const Mesh& m1, const Mesh& m2) {
    Mesh m;
    m.add(m1, m2);
    return m;
  }

}  // namespace OpenMD

#endif /*HYDRODYNAMICS_MESH_HPP_*/
