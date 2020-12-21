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

#include "applications/dynamicProps/LegendreCorrFuncZ.hpp"
#include "math/LegendrePolynomial.hpp"
#include "utils/simError.h"
#include "utils/Revision.hpp"
#include <sstream>

namespace OpenMD {
  LegendreCorrFuncZ::LegendreCorrFuncZ(SimInfo* info, 
                                       const std::string& filename, 
                                       const std::string& sele1, 
                                       const std::string& sele2, 
                                       int order, int nZbins, int axis)
    : ObjectACF<Vector3d>(info, filename, sele1, sele2,
                          DataStorage::dslPosition |
                          DataStorage::dslAmat ), nZBins_(nZbins),
    axis_(axis) {
    
    setCorrFuncType("Legendre Correlation Function of Z");
    setOutputName(getPrefix(dumpFilename_) + ".lcorrZ");
    
    std::stringstream params;
    params << " order = " << order
           << ", nzbins = " << nZBins_;
    const std::string paramString = params.str();
    setParameterString( paramString );

    if (!uniqueSelections_) {
      seleMan2_ = seleMan1_;
    }

    // Compute complementary axes to the privileged axis
    xaxis_ = (axis_ + 1) % 3;
    yaxis_ = (axis_ + 2) % 3;

    switch(axis_) {
    case 0:
      axisLabel_ = "x";
      break;
    case 1:
      axisLabel_ = "y";
      break;
    case 2:
    default:
      axisLabel_ = "z";
      break;
    }

    rotMats_.resize(nTimeBins_);
    zbin_.resize(nTimeBins_);
    histogram_.resize(nTimeBins_);
    counts_.resize(nTimeBins_);
    for (unsigned int i = 0; i < nTimeBins_; i++) {
      histogram_[i].resize(nZBins_);
      std::fill(histogram_[i].begin(), histogram_[i].end(), 0.0);
      counts_[i].resize(nZBins_);
      std::fill(counts_[i].begin(), counts_[i].end(), 0);
    }
    LegendrePolynomial polynomial(order);
    legendre_ = polynomial.getLegendrePolynomial(order);
  }

  void LegendreCorrFuncZ::computeFrame(int frame) {
    Mat3x3d hmat = currentSnapshot_ ->getHmat();
    boxZ_ = hmat(axis_,axis_);
    halfBoxZ_ = boxZ_ / 2.0;      

    ObjectACF<Vector3d>::computeFrame(frame);
  }
  
  int LegendreCorrFuncZ::computeProperty1(int frame, StuntDouble* sd) {
    
    RotMat3x3d A = sd->getA();
    rotMats_[frame].push_back( A );
    
    Vector3d pos = sd->getPos();
    if (info_->getSimParams()->getUsePeriodicBoundaryConditions())
      currentSnapshot_->wrapVector(pos);
    int zBin = int(nZBins_ * (halfBoxZ_ + pos[axis_]) / boxZ_);
    zbin_[frame].push_back(zBin);
    
    return rotMats_[frame].size() - 1;
  }

  Vector3d LegendreCorrFuncZ::calcCorrVal(int frame1, int frame2,
                                          int id1, int id2) {
    
    Vector3d v1x = rotMats_[frame1][id1].getRow(xaxis_);
    Vector3d v1y = rotMats_[frame1][id1].getRow(yaxis_);    
    Vector3d v1z = rotMats_[frame1][id1].getRow(axis_);

    Vector3d v2x = rotMats_[frame2][id2].getRow(xaxis_);
    Vector3d v2y = rotMats_[frame2][id2].getRow(yaxis_);    
    Vector3d v2z = rotMats_[frame2][id2].getRow(axis_);
    
    RealType uxprod = legendre_.evaluate(dot(v1x, v2x)/(v1x.length()*v2x.length()));
    RealType uyprod = legendre_.evaluate(dot(v1y, v2y)/(v1y.length()*v2y.length()));
    RealType uzprod = legendre_.evaluate(dot(v1z, v2z)/(v1z.length()*v2z.length()));

    return Vector3d(uxprod, uyprod, uzprod);
  }


  void LegendreCorrFuncZ::correlateFrames(int frame1, int frame2,
					  int timeBin) {
    std::vector<int> s1;
    std::vector<int> s2;

    std::vector<int>::iterator i1;
    std::vector<int>::iterator i2;

    Vector3d corrVal(0.0);

    s1 = sele1ToIndex_[frame1];

    if (uniqueSelections_)
      s2 = sele2ToIndex_[frame2];
    else
      s2 = sele1ToIndex_[frame2];

    for (i1 = s1.begin(), i2 = s2.begin();
         i1 != s1.end() && i2 != s2.end(); ++i1, ++i2){

      // If the selections are dynamic, they might not have the
      // same objects in both frames, so we need to roll either of
      // the selections until we have the same object to
      // correlate.

      while ( i1 != s1.end() && *i1 < *i2 ) {
        ++i1;
      }

      while ( i2 != s2.end() && *i2 < *i1 ) {
        ++i2;
      }

      if ( i1 == s1.end() || i2 == s2.end() ) break;

      corrVal = calcCorrVal(frame1, frame2, i1 - s1.begin(), i2 - s2.begin());
      int zBin = zbin_[frame1][ i1 - s1.begin() ];
      histogram_[timeBin][zBin] += corrVal;
      counts_[timeBin][zBin]++;

    }
  }


  void LegendreCorrFuncZ::postCorrelate() {
    for (unsigned int i =0 ; i < nTimeBins_; ++i) {
      for (unsigned int j = 0; j < nZBins_; ++j) {
        if (counts_[i][j] > 0) {
          histogram_[i][j] /= counts_[i][j];
        }
      }
    }
  }

  void LegendreCorrFuncZ::validateSelection(SelectionManager& seleMan) {
    StuntDouble* sd;
    int i;    
    for (sd = seleMan1_.beginSelected(i); sd != NULL;
         sd = seleMan1_.nextSelected(i)) {
      if (!sd->isDirectionalAtom()) {
	sprintf(painCave.errMsg,
                "LegendreCorrFunc::validateSelection Error: "
                "selected atoms are not Directional\n");
	painCave.isFatal = 1;
	simError();        
      }
    }
    
  }

  void LegendreCorrFuncZ::writeCorrelate() {
    std::ofstream ofs(getOutputFileName().c_str());

    if (ofs.is_open()) {
      Revision r;
      
      ofs << "# " << getCorrFuncType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_ ;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      ofs << "# privilegedAxis computed as " << axisLabel_ << " axis \n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";

      ofs << "#time\tPn(costheta_z)\n";

      for (unsigned int i = 0; i < nTimeBins_; ++i) {

        ofs << times_[i]-times_[0];

        for (unsigned int j = 0; j < nZBins_; ++j) {          
          ofs << "\t" << histogram_[i][j](2);
        }
        ofs << "\n";
      }
            
    } else {
      sprintf(painCave.errMsg,
              "LegendreCorrFuncZ::writeCorrelate Error: failed to open %s\n",
              getOutputFileName().c_str());
      painCave.isFatal = 1;
      simError();        
    }
    ofs.close();    
  }
}
