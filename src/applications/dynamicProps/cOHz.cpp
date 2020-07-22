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

#include "applications/dynamicProps/cOHz.hpp"
#include "math/LegendrePolynomial.hpp"
#include "utils/simError.h"
#include "utils/Revision.hpp"
#include <sstream>

namespace OpenMD {
  COHZ::COHZ(SimInfo* info, 
             const std::string& filename, 
             const std::string& sele1, 
             const std::string& sele2, 
             int order, int nZbins, int axis)
    : MoleculeACF<Vector<RealType, 4> >(info, filename, sele1, sele2,
                                        DataStorage::dslPosition |
                                        DataStorage::dslAmat),
    nZBins_(nZbins), axis_(axis) {

    setCorrFuncType("Legendre Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".lcorr");

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
      std::fill(histogram_[i].begin(), histogram_[i].end(),
                Vector<RealType, 4>(0.0) );
      counts_[i].resize(nZBins_);
      std::fill(counts_[i].begin(), counts_[i].end(), 0);
    }
    LegendrePolynomial polynomial(order);
    legendre_ = polynomial.getLegendrePolynomial(order);
  }

  void COHZ::computeFrame(int frame) {
    Mat3x3d hmat = currentSnapshot_ ->getHmat();
    boxZ_ = hmat(axis_,axis_);
    halfBoxZ_ = boxZ_ / 2.0;      

    MoleculeACF<Vector<RealType, 4> >::computeFrame(frame);
  }

  int COHZ::computeProperty1(int frame, Molecule* mol) {
    
    RotMat3x3d A = mol->getRigidBodyAt(0)->getA();
    rotMats_[frame].push_back( A );
    
    Vector3d pos = mol->getCom();
    if (info_->getSimParams()->getUsePeriodicBoundaryConditions())
      currentSnapshot_->wrapVector(pos);
    int zBin = int(nZBins_ * (halfBoxZ_ + pos[axis_]) / boxZ_);
    zbin_[frame].push_back(zBin);
    
    return rotMats_[frame].size() - 1;
  }
  
  Vector<RealType, 4> COHZ::calcCorrVal(int frame1, int frame2,
                                        int id1, int id2) {
    
    // Vectors v1x, v1y, and v1z are the body-fixed axes on the
    // molecule in frame 1 in the laboratory frame.

    // Vectors v2x, v2y, and v2z are the body-fixed axes on the
    // molecule in frame 2 in the laboratory frame.

    // Vectors u1 & u2 are the first OH bond vector in frames 1 & 2
    // respectively.  Here we assume SPC/E geometry.

    // Vectors w1 & w2 are the second OH bond vector in frames 1 & 2
    // respectively.  Here we assume SPC/E geometry.

    // Vectors h1 & h2 are the HH bond vectors in frames 1 & 2
    // respectively.  Here we assume SPC/E geometry again.

    // Vector3d v1x = sd1->getA(frame1).getRow(0);
    // Vector3d v2x = sd2->getA(frame2).getRow(0);

    Vector3d v1y = rotMats_[frame1][id1].getRow(yaxis_);    
    Vector3d v1z = rotMats_[frame1][id1].getRow(axis_);

    Vector3d v2y = rotMats_[frame2][id2].getRow(yaxis_);    
    Vector3d v2z = rotMats_[frame2][id2].getRow(axis_);

    Vector3d u1 = 0.81649 * v1y + 0.57736 * v1z;
    Vector3d u2 = 0.81649 * v2y + 0.57736 * v2z;

    Vector3d w1 = -0.81649 * v1y + 0.57736 * v1z;
    Vector3d w2 = -0.81649 * v2y + 0.57736 * v2z;

    Vector3d h1 = 1.63298 * v1y;
    Vector3d h2 = 1.63298 * v2y;

    // result is a Vector<RealType, 4> with Dipole, OH1, OH2, and HH:
    Vector<RealType, 4> r(0.0);
    r[0] = legendre_.evaluate(dot(v1z, v2z)/(v1z.length()*v2z.length()));
    r[1] = legendre_.evaluate(dot(u1, u2)/(u1.length()*u2.length()));
    r[2] = legendre_.evaluate(dot(w1, w2)/(w1.length()*w2.length()));
    r[3] = legendre_.evaluate(dot(h1, h2)/(h1.length()*h2.length()));
    return r;
  }

  void COHZ::correlateFrames(int frame1, int frame2,
                             int timeBin) {
    std::vector<int> s1;
    std::vector<int> s2;
    
    std::vector<int>::iterator i1;
    std::vector<int>::iterator i2;
    
    Vector<RealType, 4> corrVal(0.0);
    
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
  
  void COHZ::postCorrelate() {
    for (unsigned int i =0 ; i < nTimeBins_; ++i) {
      for (unsigned int j = 0; j < nZBins_; ++j) {
        if (counts_[i][j] > 0) {
          histogram_[i][j] /= counts_[i][j];          
        }
      }
    }
  }

  void COHZ::validateSelection(SelectionManager& seleMan) {
    Molecule* mol;
    int i;    
    for (mol = seleMan1_.beginSelectedMolecule(i); mol != NULL;
         mol = seleMan1_.nextSelectedMolecule(i)) {
      if (mol->getNRigidBodies() < 1) {
        sprintf(painCave.errMsg,
                "COHZ::validateSelection Error: "
                "at least one selected molecule does not have a rigid body\n");
        painCave.isFatal = 1;
        simError();        
      }
    }    
  }

  void COHZ::writeCorrelate() {
    
    std::string Dfile = getOutputFileName() + "D";
    std::string OHfile = getOutputFileName() + "OH";
    std::string HHfile = getOutputFileName() + "HH";
    
    std::ofstream ofs1(Dfile.c_str());
    std::ofstream ofs2(OHfile.c_str());
    std::ofstream ofs3(HHfile.c_str());

    if (ofs1.is_open()) {
      Revision r;
      
      ofs1 << "# " << getCorrFuncType() << " for dipole vectors in water\n";
      ofs1 << "# OpenMD " << r.getFullRevision() << "\n";
      ofs1 << "# " << r.getBuildDate() << "\n";
      ofs1 << "# selection script1: \"" << selectionScript1_ ;
      ofs1 << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      ofs1 << "# privilegedAxis computed as " << axisLabel_ << " axis \n";
      if (!paramString_.empty())
        ofs1 << "# parameters: " << paramString_ << "\n";

      ofs1 << "#time\tPn(costheta_z)\n";

      for (unsigned int i = 0; i < nTimeBins_; ++i) {

        ofs1 << times_[i]-times_[0];

        for (unsigned int j = 0; j < nZBins_; ++j) {          
          ofs1 << "\t" << histogram_[i][j][0];
        }
        ofs1 << "\n";
      }
            
    } else {
      sprintf(painCave.errMsg,
              "COHz::writeCorrelate Error: failed to open %s\n",
              Dfile.c_str());
      painCave.isFatal = 1;
      simError();        
    }
    ofs1.close();    
  
    if (ofs2.is_open()) {
      Revision r;
      
      ofs2 << "# " << getCorrFuncType() << " for OH bond vectors in water\n";
      ofs2 << "# OpenMD " << r.getFullRevision() << "\n";
      ofs2 << "# " << r.getBuildDate() << "\n";
      ofs2 << "# selection script1: \"" << selectionScript1_ ;
      ofs2 << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      ofs2 << "# privilegedAxis computed as " << axisLabel_ << " axis \n";
      if (!paramString_.empty())
        ofs2 << "# parameters: " << paramString_ << "\n";
      
      ofs2 << "#time\tPn(costheta_z)\n";
      
      for (unsigned int i = 0; i < nTimeBins_; ++i) {
        
        ofs2 << times_[i]-times_[0];
        
        for (unsigned int j = 0; j < nZBins_; ++j) {          
          ofs2 << "\t" << 0.5 * (histogram_[i][j][1] + histogram_[i][j][2]);
        }
        ofs2 << "\n";
      }
      
    } else {
      sprintf(painCave.errMsg,
              "COHz::writeCorrelate Error: failed to open %s\n",
              OHfile.c_str());
      painCave.isFatal = 1;
      simError();        
    }
    ofs2.close();

    if (ofs3.is_open()) {
      Revision r;
      
      ofs3 << "# " << getCorrFuncType() << " for HH bond vectors in water\n";
      ofs3 << "# OpenMD " << r.getFullRevision() << "\n";
      ofs3 << "# " << r.getBuildDate() << "\n";
      ofs3 << "# selection script1: \"" << selectionScript1_ ;
      ofs3 << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      ofs3 << "# privilegedAxis computed as " << axisLabel_ << " axis \n";
      if (!paramString_.empty())
        ofs3 << "# parameters: " << paramString_ << "\n";
      
      ofs3 << "#time\tPn(costheta_z)\n";
      
      for (unsigned int i = 0; i < nTimeBins_; ++i) {
        
        ofs3 << times_[i]-times_[0];
        
        for (unsigned int j = 0; j < nZBins_; ++j) {          
          ofs3 << "\t" << histogram_[i][j][3];
        }
        ofs3 << "\n";
      }
      
    } else {
      sprintf(painCave.errMsg,
              "COHz::writeCorrelate Error: failed to open %s\n",
              HHfile.c_str());
      painCave.isFatal = 1;
      simError();        
    }
    ofs3.close();    
  }
}
