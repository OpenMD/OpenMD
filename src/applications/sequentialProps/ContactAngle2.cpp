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

#include "applications/sequentialProps/ContactAngle2.hpp"

#include <algorithm>
#include <functional>
#include <sstream>

#include "io/DumpReader.hpp"
#include "math/Eigenvalue.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  ContactAngle2::ContactAngle2(SimInfo* info, const std::string& filename,
                               const std::string& sele1,
                               const std::string& sele2, RealType solidZ,
                               RealType centroidX, RealType centroidY,
                               RealType threshDens, RealType bufferLength,
                               int nrbins, int nzbins) :
      SequentialAnalyzer(info, filename, sele1, sele2),
      solidZ_(solidZ), centroidX_(centroidX), centroidY_(centroidY),
      threshDens_(threshDens), bufferLength_(bufferLength), nRBins_(nrbins),
      nZBins_(nzbins) {
    setOutputName(getPrefix(filename) + ".ca2");

    std::stringstream params;
    params << " referenceZ = " << solidZ_ << ", centroid = (" << centroidX_
           << ", " << centroidY_ << ")"
           << ", threshDens = " << threshDens_
           << ", bufferLength = " << bufferLength_ << ", nbins = " << nRBins_
           << ", nbins_z = " << nZBins_;

    const std::string paramString = params.str();
    setParameterString(paramString);
  }

  void ContactAngle2::doFrame(int) {
    StuntDouble* sd;
    int i;

    // set up the bins for density analysis

    Mat3x3d hmat = info_->getSnapshotManager()->getCurrentSnapshot()->getHmat();
    RealType len = std::min(hmat(0, 0), hmat(1, 1));
    RealType zLen = hmat(2, 2);

    RealType dr = len / (RealType)nRBins_;
    RealType dz = zLen / (RealType)nZBins_;

    std::vector<std::vector<RealType>> histo;
    histo.resize(nRBins_);
    for (unsigned int i = 0; i < histo.size(); ++i) {
      histo[i].resize(nZBins_);
      std::fill(histo[i].begin(), histo[i].end(), 0.0);
    }

    if (evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    Vector3d com(centroidX_, centroidY_, solidZ_);

    // now that we have the centroid, we can make cylindrical density maps
    Vector3d pos;
    RealType r;
    RealType z;

    for (sd = seleMan1_.beginSelected(i); sd != NULL;
         sd = seleMan1_.nextSelected(i)) {
      pos = sd->getPos() - com;

      // r goes from zero upwards
      r = sqrt(pow(pos.x(), 2) + pow(pos.y(), 2));
      // z is possibly symmetric around 0
      z = pos.z();

      int whichRBin = int(r / dr);
      int whichZBin = int((zLen / 2.0 + z) / dz);

      if ((whichRBin < int(nRBins_)) && (whichZBin >= 0) &&
          (whichZBin < int(nZBins_))) {
        histo[whichRBin][whichZBin] += sd->getMass();
      }
    }

    for (unsigned int i = 0; i < histo.size(); ++i) {
      RealType rL       = i * dr;
      RealType rU       = rL + dr;
      RealType volSlice = Constants::PI * dz * ((rU * rU) - (rL * rL));

      for (unsigned int j = 0; j < histo[i].size(); ++j) {
        histo[i][j] *= Constants::densityConvert / volSlice;
      }
    }

    std::vector<Vector<RealType, 2>> points;
    points.clear();

    for (unsigned int j = 0; j < nZBins_; ++j) {
      // The z coordinates were measured relative to the selection
      // center of mass.  However, we're interested in the elevation
      // above the solid surface.  Also, the binning was done around
      // zero with enough bins to cover the zLength of the box:

      RealType thez    = com.z() - solidZ_ - zLen / 2.0 + dz * (j + 0.5);
      bool aboveThresh = false;
      bool foundThresh = false;
      int rloc         = 0;

      for (std::size_t i = 0; i < nRBins_; ++i) {
        if (histo[i][j] >= threshDens_) aboveThresh = true;

        if (aboveThresh && (histo[i][j] <= threshDens_)) {
          rloc        = i;
          foundThresh = true;
          aboveThresh = false;
        }
      }
      if (foundThresh) {
        Vector<RealType, 2> point;
        point[0] = dr * (rloc + 0.5);
        point[1] = thez;

        if (thez > bufferLength_) { points.push_back(point); }
      }
    }

    int numPoints = points.size();

    // Compute the average of the data points.
    Vector<RealType, 2> average = points[0];
    int i0;
    for (i0 = 1; i0 < numPoints; ++i0) {
      average += points[i0];
    }
    RealType invNumPoints = ((RealType)1) / (RealType)numPoints;
    average *= invNumPoints;

    DynamicRectMatrix<RealType> mat(4, 4);
    int row, col;
    for (row = 0; row < 4; ++row) {
      for (col = 0; col < 4; ++col) {
        mat(row, col) = 0.0;
      }
    }
    for (int i = 0; i < numPoints; ++i) {
      RealType x   = points[i][0];
      RealType y   = points[i][1];
      RealType x2  = x * x;
      RealType y2  = y * y;
      RealType xy  = x * y;
      RealType r2  = x2 + y2;
      RealType xr2 = x * r2;
      RealType yr2 = y * r2;
      RealType r4  = r2 * r2;

      mat(0, 1) += x;
      mat(0, 2) += y;
      mat(0, 3) += r2;
      mat(1, 1) += x2;
      mat(1, 2) += xy;
      mat(1, 3) += xr2;
      mat(2, 2) += y2;
      mat(2, 3) += yr2;
      mat(3, 3) += r4;
    }
    mat(0, 0) = (RealType)numPoints;

    for (row = 0; row < 4; ++row) {
      for (col = 0; col < row; ++col) {
        mat(row, col) = mat(col, row);
      }
    }

    for (row = 0; row < 4; ++row) {
      for (col = 0; col < 4; ++col) {
        mat(row, col) *= invNumPoints;
      }
    }

    JAMA::Eigenvalue<RealType> eigensystem(mat);
    DynamicRectMatrix<RealType> evects(4, 4);
    DynamicVector<RealType> evals(4);

    eigensystem.getRealEigenvalues(evals);
    eigensystem.getV(evects);

    DynamicVector<RealType> evector = evects.getColumn(0);
    RealType inv = ((RealType)1) / evector[3];  // beware zero divide
    RealType coeff[3];
    for (row = 0; row < 3; ++row) {
      coeff[row] = inv * evector[row];
    }

    Vector<RealType, 2> center;

    center[0] = -((RealType)0.5) * coeff[1];
    center[1] = -((RealType)0.5) * coeff[2];
    RealType radius =
        sqrt(fabs(center[0] * center[0] + center[1] * center[1] - coeff[0]));

    int i1;
    for (i1 = 0; i1 < 100; ++i1) {
      // Update the iterates.
      Vector<RealType, 2> current = center;

      // Compute average L, dL/da, dL/db.
      RealType lenAverage               = (RealType)0;
      Vector<RealType, 2> derLenAverage = Vector<RealType, 2>(0.0);
      for (i0 = 0; i0 < numPoints; ++i0) {
        Vector<RealType, 2> diff = points[i0] - center;
        RealType length          = diff.length();
        if (length > 1e-6) {
          lenAverage += length;
          RealType invLength = ((RealType)1) / length;
          derLenAverage -= invLength * diff;
        }
      }
      lenAverage *= invNumPoints;
      derLenAverage *= invNumPoints;

      center = average + lenAverage * derLenAverage;
      radius = lenAverage;

      Vector<RealType, 2> diff = center - current;
      if (fabs(diff[0]) <= 1e-6 && fabs(diff[1]) <= 1e-6) { break; }
    }

    RealType zCen  = center[1];
    RealType rDrop = radius;
    RealType ca;

    if (fabs(zCen) > rDrop) {
      ca = 180.0;
    } else {
      ca = 90.0 + asin(zCen / rDrop) * (180.0 / Constants::PI);
    }

    values_.push_back(ca);
  }
}  // namespace OpenMD
