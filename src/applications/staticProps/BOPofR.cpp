/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

#include "applications/staticProps/BOPofR.hpp"

#include "brains/Thermo.hpp"
#include "io/DumpReader.hpp"
#include "math/Wigner3jm.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/Revision.hpp"
#include "utils/simError.h"

using namespace MATPACK;
namespace OpenMD {

  BOPofR::BOPofR(SimInfo* info, const std::string& filename,
                 const std::string& sele, double rCut, unsigned int nbins,
                 RealType len) :
      StaticAnalyser(info, filename, nbins),
      selectionScript_(sele), seleMan_(info), evaluator_(info) {
    setOutputName(getPrefix(filename) + ".bo");
    setAnalysisType("Bond Order Parameter(r)");

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // Set up cutoff radius and order of the Legendre Polynomial:

    rCut_ = rCut;
    len_  = len;

    std::stringstream params;
    params << " rcut = " << rCut_ << ", len = " << len_
           << ", nbins = " << nBins_;
    const std::string paramString = params.str();
    setParameterString(paramString);

    deltaR_ = len_ / nBins_;
    RCount_.resize(nBins_);
    WofR_.resize(nBins_);
    QofR_.resize(nBins_);

    for (unsigned int i = 0; i < nBins_; i++) {
      RCount_[i] = 0;
      WofR_[i]   = 0;
      QofR_[i]   = 0;
    }

    // Make arrays for Wigner3jm
    RealType* THRCOF = new RealType[2 * lMax_ + 1];
    // Variables for Wigner routine
    RealType lPass, m1Pass, m2m, m2M;
    int error, mSize;
    mSize = 2 * lMax_ + 1;

    for (int l = 0; l <= lMax_; l++) {
      lPass = (RealType)l;
      for (int m1 = -l; m1 <= l; m1++) {
        m1Pass = (RealType)m1;

        std::pair<int, int> lm = std::make_pair(l, m1);

        // Zero work array
        for (int ii = 0; ii < 2 * l + 1; ii++) {
          THRCOF[ii] = 0.0;
        }

        // Get Wigner coefficients
        Wigner3jm(lPass, lPass, lPass, m1Pass, m2m, m2M, THRCOF, mSize, error);

        m2Min[lm] = (int)floor(m2m);
        m2Max[lm] = (int)floor(m2M);

        for (int mmm = 0; mmm <= (int)(m2M - m2m); mmm++) {
          w3j[lm].push_back(THRCOF[mmm]);
        }
      }
    }

    delete[] THRCOF;
    THRCOF = NULL;
  }

  void BOPofR::initializeHistogram() {
    for (unsigned int i = 0; i < nBins_; i++) {
      RCount_[i] = 0;
      WofR_[i]   = 0;
      QofR_[i]   = 0;
    }
  }

  void BOPofR::process() {
    Molecule* mol;
    Atom* atom;
    int myIndex;
    SimInfo::MoleculeIterator mi;
    Molecule::AtomIterator ai;
    StuntDouble* sd;
    Vector3d vec;
    RealType costheta;
    RealType phi;
    RealType r;
    Vector3d rCOM;
    RealType distCOM;
    Vector3d pos;
    Vector3d CenterOfMass;
    std::map<std::pair<int, int>, ComplexType> q;
    std::vector<RealType> q_l;
    std::vector<RealType> q2;
    std::vector<ComplexType> w;
    std::vector<ComplexType> w_hat;
    std::vector<RealType> Q2;
    std::vector<RealType> Q;
    std::vector<ComplexType> W;
    std::vector<ComplexType> W_hat;
    int nBonds;
    SphericalHarmonic sphericalHarmonic;
    int i;
    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);
    int nFrames   = reader.getNFrames();
    frameCounter_ = 0;

    Thermo thermo(info_);

    q_l.resize(lMax_ + 1);
    q2.resize(lMax_ + 1);
    w.resize(lMax_ + 1);
    w_hat.resize(lMax_ + 1);

    Q2.resize(lMax_ + 1);
    Q.resize(lMax_ + 1);
    W.resize(lMax_ + 1);
    W_hat.resize(lMax_ + 1);

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      frameCounter_++;
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      CenterOfMass     = thermo.getCom();
      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      // outer loop is over the selected StuntDoubles:

      for (sd = seleMan_.beginSelected(i); sd != NULL;
           sd = seleMan_.nextSelected(i)) {
        myIndex = sd->getGlobalIndex();

        nBonds = 0;

        for (int l = 0; l <= lMax_; l++) {
          for (int m = -l; m <= l; m++) {
            q[std::make_pair(l, m)] = 0.0;
          }
        }
        pos  = sd->getPos();
        rCOM = CenterOfMass - pos;
        if (usePeriodicBoundaryConditions_) currentSnapshot_->wrapVector(rCOM);
        distCOM = rCOM.length();

        // inner loop is over all other atoms in the system:

        for (mol = info_->beginMolecule(mi); mol != NULL;
             mol = info_->nextMolecule(mi)) {
          for (atom = mol->beginAtom(ai); atom != NULL;
               atom = mol->nextAtom(ai)) {
            if (atom->getGlobalIndex() != myIndex) {
              vec = pos - atom->getPos();

              if (usePeriodicBoundaryConditions_)
                currentSnapshot_->wrapVector(vec);

              // Calculate "bonds" and build Q_lm(r) where
              //      Q_lm = Y_lm(theta(r),phi(r))
              // The spherical harmonics are wrt any arbitrary coordinate
              // system, we choose standard spherical coordinates

              r = vec.length();

              // Check to see if neighbor is in bond cutoff

              if (r < rCut_) {
                costheta = vec.z() / r;
                phi      = atan2(vec.y(), vec.x());

                for (int l = 0; l <= lMax_; l++) {
                  sphericalHarmonic.setL(l);
                  for (int m = -l; m <= l; m++) {
                    sphericalHarmonic.setM(m);
                    q[std::make_pair(l, m)] +=
                        sphericalHarmonic.getValueAt(costheta, phi);
                  }
                }
                nBonds++;
              }
            }
          }
        }

        for (int l = 0; l <= lMax_; l++) {
          q2[l] = 0.0;
          for (int m = -l; m <= l; m++) {
            q[std::make_pair(l, m)] /= (RealType)nBonds;
            q2[l] += norm(q[std::make_pair(l, m)]);
          }
          q_l[l] = sqrt(q2[l] * 4.0 * Constants::PI / (RealType)(2 * l + 1));
        }

        // Find Third Order Invariant W_l

        for (int l = 0; l <= lMax_; l++) {
          w[l] = 0.0;
          for (int m1 = -l; m1 <= l; m1++) {
            std::pair<int, int> lm = std::make_pair(l, m1);
            for (int mmm = 0; mmm <= (m2Max[lm] - m2Min[lm]); mmm++) {
              int m2 = m2Min[lm] + mmm;
              int m3 = -m1 - m2;
              w[l] += w3j[lm][mmm] * q[lm] * q[std::make_pair(l, m2)] *
                      q[std::make_pair(l, m3)];
            }
          }

          w_hat[l] = w[l] / pow(q2[l], RealType(1.5));
        }

        collectHistogram(q_l, w_hat, distCOM);

        //  printf( "%s  %18.10g %18.10g %18.10g %18.10g \n",
        //  sd->getType().c_str(),pos[0],pos[1],pos[2],real(w_hat[6]));
      }
    }

    writeOrderParameter();
  }

  IcosahedralOfR::IcosahedralOfR(SimInfo* info, const std::string& filename,
                                 const std::string& sele, double rCut,
                                 unsigned int nbins, RealType len) :
      BOPofR(info, filename, sele, rCut, nbins, len) {
    setAnalysisType("Icosahedral Bond Order Parameter(r)");
  }

  void IcosahedralOfR::collectHistogram(std::vector<RealType> q,
                                        std::vector<ComplexType> what,
                                        RealType distCOM) {
    if (distCOM < len_) {
      // Figure out where this distance goes...
      int whichBin = int(distCOM / deltaR_);
      RCount_[whichBin]++;

      if (real(what[6]) < -0.15) { WofR_[whichBin]++; }
      if (q[6] > 0.5) { QofR_[whichBin]++; }
    }
  }

  FCCOfR::FCCOfR(SimInfo* info, const std::string& filename,
                 const std::string& sele, double rCut, unsigned int nbins,
                 RealType len) :
      BOPofR(info, filename, sele, rCut, nbins, len) {
    setAnalysisType("FCC Bond Order Parameter(r)");
  }

  void FCCOfR::collectHistogram(std::vector<RealType>,
                                std::vector<ComplexType> what,
                                RealType distCOM) {
    if (distCOM < len_) {
      // Figure out where this distance goes...
      int whichBin = int(distCOM / deltaR_);
      RCount_[whichBin]++;

      if (real(what[4]) < -0.12) { WofR_[whichBin]++; }
    }
  }

  void IcosahedralOfR::writeOrderParameter() {
    Revision rev;
    std::ofstream osq((getOutputFileName() + "qr").c_str());

    if (osq.is_open()) {
      osq << "# " << getAnalysisType() << "\n";
      osq << "# OpenMD " << rev.getFullRevision() << "\n";
      osq << "# " << rev.getBuildDate() << "\n";
      osq << "# selection script: \"" << selectionScript_ << "\"\n";
      if (!paramString_.empty())
        osq << "# parameters: " << paramString_ << "\n";

      // Normalize by number of frames and write it out:

      for (unsigned int i = 0; i < nBins_; ++i) {
        RealType Rval = (i + 0.5) * deltaR_;
        osq << Rval;
        if (RCount_[i] == 0) {
          osq << "\t" << 0;
          osq << "\n";
        } else {
          osq << "\t" << (RealType)QofR_[i] / (RealType)RCount_[i];
          osq << "\n";
        }
      }

      osq.close();

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "IcosahedralOfR: unable to open %s\n",
               (getOutputFileName() + "q").c_str());
      painCave.isFatal = 1;
      simError();
    }

    std::ofstream osw((getOutputFileName() + "wr").c_str());

    if (osw.is_open()) {
      osw << "# " << getAnalysisType() << "\n";
      osw << "# OpenMD " << rev.getFullRevision() << "\n";
      osw << "# " << rev.getBuildDate() << "\n";
      osw << "# selection script: \"" << selectionScript_ << "\"\n";
      if (!paramString_.empty())
        osw << "# parameters: " << paramString_ << "\n";

      // Normalize by number of frames and write it out:
      for (unsigned int i = 0; i < nBins_; ++i) {
        RealType Rval = deltaR_ * (i + 0.5);
        osw << Rval;
        if (RCount_[i] == 0) {
          osw << "\t" << 0;
          osw << "\n";
        } else {
          osw << "\t" << (RealType)WofR_[i] / (RealType)RCount_[i];
          osw << "\n";
        }
      }

      osw.close();
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "IcosahedralOfR: unable to open %s\n",
               (getOutputFileName() + "w").c_str());
      painCave.isFatal = 1;
      simError();
    }
  }
  void FCCOfR::writeOrderParameter() {
    std::ofstream osw((getOutputFileName() + "wr").c_str());

    if (osw.is_open()) {
      Revision rev;
      osw << "# " << getAnalysisType() << "\n";
      osw << "# OpenMD " << rev.getFullRevision() << "\n";
      osw << "# " << rev.getBuildDate() << "\n";
      osw << "# selection script: \"" << selectionScript_ << "\"\n";
      if (!paramString_.empty())
        osw << "# parameters: " << paramString_ << "\n";

      // Normalize by number of frames and write it out:
      for (unsigned int i = 0; i < nBins_; ++i) {
        RealType Rval = deltaR_ * (i + 0.5);
        osw << Rval;
        if (RCount_[i] == 0) {
          osw << "\t" << 0;
          osw << "\n";
        } else {
          osw << "\t" << (RealType)WofR_[i] / (RealType)RCount_[i];
          osw << "\n";
        }
      }

      osw.close();
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "FCCOfR: unable to open %s\n",
               (getOutputFileName() + "w").c_str());
      painCave.isFatal = 1;
      simError();
    }
  }
}  // namespace OpenMD
