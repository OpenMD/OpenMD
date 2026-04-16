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

#include "applications/sequentialProps/Qsurf.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/Revision.hpp"
#include "utils/simError.h"

namespace OpenMD {

  Qsurf::Qsurf(SimInfo* info, const std::string& filename,
               const std::string& sele1, const std::string& sele2,
               RealType rCut, RealType voxelSize, RealType gaussWidth) :
    SequentialAnalyzer(info, filename, sele1, sele2),
    rCut_(rCut), voxelSize_(voxelSize), gaussWidth_(gaussWidth),
    nBinsX_(0), nBinsY_(0), nBinsZ_(0) {
    setSequenceType("Qsurf");

    std::ostringstream params;
    params << "rCut = " << rCut_ << ", voxelSize = " << voxelSize_
           << ", gaussWidth = " << gaussWidth_;
    setParameterString(params.str());

    // Determine grid dimensions from the initial snapshot's box
    Mat3x3d hmat =
      info->getSnapshotManager()->getCurrentSnapshot()->getHmat();
    nBinsX_ = int(hmat(0, 0) / voxelSize_);
    nBinsY_ = int(hmat(1, 1) / voxelSize_);
    nBinsZ_ = int(hmat(2, 2) / voxelSize_);
  }

  void Qsurf::doSequence() {
    preSequence();

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();

    bool usePeriodicBoundaryConditions =
      info_->getSimParams()->getUsePeriodicBoundaryConditions();

    // Gaussian kernel truncation radius (in voxel units)
    int kMax   = int(5.0 * gaussWidth_ / voxelSize_);
    int kSqLim = kMax * kMax;

    // Determine the base name for output files.
    // If outputFilename_ is "foo.dx", baseName is "foo"
    // If outputFilename_ is "foo", baseName is "foo"
    std::string baseName = outputFilename_;
    std::string::size_type pos = baseName.rfind(".dx");
    if (pos != std::string::npos && pos == baseName.length() - 3) {
      baseName = baseName.substr(0, pos);
    }

    int nFramesWritten = 0;

    for (frame_ = 0; frame_ < nFrames; frame_ += step_) {
      reader.readFrame(frame_);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      RealType time     = currentSnapshot_->getTime();
      Mat3x3d hmat      = currentSnapshot_->getHmat();
      Vector3d halfBox =
	Vector3d(hmat(0, 0), hmat(1, 1), hmat(2, 2)) / 2.0;

      if (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
      }
      if (evaluator2_.isDynamic()) {
        seleMan2_.setSelectionSet(evaluator2_.evaluate());
      }

      // ----------------------------------------------------------
      // 1) Allocate and zero the 3D grids for this frame
      // ----------------------------------------------------------
      std::vector<std::vector<std::vector<RealType>>> qField(nBinsX_, std::vector<std::vector<RealType>>(nBinsY_, std::vector<RealType>(nBinsZ_, 0.0)));

      std::vector<std::vector<std::vector<RealType>>> weight(nBinsX_, std::vector<std::vector<RealType>>(nBinsY_, std::vector<RealType>(nBinsZ_, 0.0)));

      // ----------------------------------------------------------
      // 2) Compute Q_k for each molecule in sele1, smear onto grid
      // ----------------------------------------------------------
      int isd1, isd2;
      StuntDouble* sd;
      StuntDouble* sd2;

      for (sd = seleMan1_.beginSelected(isd1); sd != NULL;
           sd = seleMan1_.nextSelected(isd1)) {
        int myIndex = sd->getGlobalIndex();

        // Find neighbors from sele2 within rCut
        std::vector<std::pair<RealType, StuntDouble*>> neighbors;

        for (sd2 = seleMan2_.beginSelected(isd2); sd2 != NULL;
             sd2 = seleMan2_.nextSelected(isd2)) {
          if (sd2->getGlobalIndex() == myIndex) continue;

          Vector3d vec = sd->getPos() - sd2->getPos();
          if (usePeriodicBoundaryConditions)
            currentSnapshot_->wrapVector(vec);

          RealType r = vec.length();
          if (r < rCut_) { neighbors.push_back(std::make_pair(r, sd2)); }
        }

        std::sort(neighbors.begin(), neighbors.end());
        int nbors = std::min((int)neighbors.size(), 4);
        int nang  = nbors * (nbors - 1) / 2;
        if (nang == 0) continue;

        // Compute Q_k (Errington-Debenedetti rescaled tetrahedrality)
        Vector3d rk = sd->getPos();
        RealType Qk = 1.0;

        for (int i = 0; i < nbors - 1; i++) {
          Vector3d rik = rk - neighbors[i].second->getPos();
          if (usePeriodicBoundaryConditions)
            currentSnapshot_->wrapVector(rik);
          rik.normalize();

          for (int j = i + 1; j < nbors; j++) {
            Vector3d rkj = rk - neighbors[j].second->getPos();
            if (usePeriodicBoundaryConditions)
              currentSnapshot_->wrapVector(rkj);
            rkj.normalize();

            RealType cospsi = dot(rik, rkj);
            Qk -= (pow(cospsi + 1.0 / 3.0, 2) * 2.25 / nang);
          }
        }

        // Smear Q_k onto the voxel grid with a Gaussian kernel
        if (usePeriodicBoundaryConditions)
          currentSnapshot_->wrapVector(rk);

        Vector3d gridPos = rk + halfBox;
        Vector3i whichVoxel(int(gridPos[0] / voxelSize_),
                            int(gridPos[1] / voxelSize_),
                            int(gridPos[2] / voxelSize_));

        RealType denom = pow(2.0 * sqrt(Constants::PI) * gaussWidth_, 3);

        for (int l = -kMax; l <= kMax; l++) {
          for (int m = -kMax; m <= kMax; m++) {
            for (int n = -kMax; n <= kMax; n++) {
              if (l * l + m * m + n * n > kSqLim) continue;

              int ll = (whichVoxel[0] + l) % nBinsX_;
              if (ll < 0) ll += nBinsX_;
              int mm = (whichVoxel[1] + m) % nBinsY_;
              if (mm < 0) mm += nBinsY_;
              int nn = (whichVoxel[2] + n) % nBinsZ_;
              if (nn < 0) nn += nBinsZ_;

              Vector3d bPos = Vector3d(ll, mm, nn) * voxelSize_ - halfBox;
              Vector3d d    = bPos - rk;
              currentSnapshot_->wrapVector(d);

              RealType exponent = -dot(d, d) / pow(2.0 * gaussWidth_, 2);
              RealType w        = exp(exponent) / denom;

              weight[ll][mm][nn] += w;
              qField[ll][mm][nn] += w * Qk;
            }
          }
        }
      }  // end sele1 loop

      // ----------------------------------------------------------
      // 3) Normalize: qField[i][j][k] = weighted average Q
      // ----------------------------------------------------------
      for (int i = 0; i < nBinsX_; i++) {
        for (int j = 0; j < nBinsY_; j++) {
          for (int k = 0; k < nBinsZ_; k++) {
            if (weight[i][j][k] > 0.0)
              qField[i][j][k] /= weight[i][j][k];
          }
        }
      }

      // ----------------------------------------------------------
      // 4) Write this frame's Q field as an OpenDX file
      // ----------------------------------------------------------
      std::ostringstream dxFileName;
      dxFileName << baseName << "." << std::setfill('0') << std::setw(4)
                 << nFramesWritten << ".dx";

      Vector3d origin = -halfBox;
      writeDXFrame(dxFileName.str(), qField, origin, time);
      nFramesWritten++;

    }  // end frame loop

    // ----------------------------------------------------------
    // 5) Write a VMD Tcl loader script
    // ----------------------------------------------------------
    std::string scriptName = baseName + ".vmd";
    writeVMDScript(scriptName, baseName, nFramesWritten);
  }

  void Qsurf::writeDXFrame(
			   const std::string& fname,
			   const std::vector<std::vector<std::vector<RealType>>>& field,
			   const Vector3d& origin,
			   RealType time) {

    std::ofstream dxStream(fname.c_str());
    if (!dxStream.is_open()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Qsurf: unable to open %s for writing\n", fname.c_str());
      painCave.isFatal = 1;
      simError();
    }

    int nItems = nBinsX_ * nBinsY_ * nBinsZ_;

    // ----------------------------------------------------------
    // OpenDX header (APBS-compatible, VMD-readable)
    //
    // To load this file in VMD:
    //   mol addfile myfile.dx type dx waitfor all
    //   mol representation Isosurface 0.85 0 0 0 1 1
    //   mol addrep top
    //
    // To animate a series, use the companion .vmd script:
    //   vmd -e basename.vmd
    //
    // In Jmol, use:
    //   isosurface cutoff 0.85 "myfile.dx"
    // ----------------------------------------------------------

    Revision r;
    dxStream << "# OpenDX Data Explorer format\n";
    dxStream << "# Generated by OpenMD Qsurf (SequentialProps)\n";
    dxStream << "# OpenMD " << r.getFullRevision() << "\n";
    dxStream << "# " << r.getBuildDate() << "\n";
    dxStream << "# Tetrahedral order parameter Q field\n";
    dxStream << "# Simulation time: " << time << " fs\n";
    dxStream << "# Grid: " << nBinsX_ << " x " << nBinsY_ << " x "
             << nBinsZ_ << "  voxelSize = " << voxelSize_ << " Ang\n";
    dxStream << "#\n";
    dxStream << "# To load in VMD:\n";
    dxStream << "#   mol new " << fname << " type dx waitfor all\n";
    dxStream << "#   mol representation Isosurface 0.85 0 0 0 1 1\n";
    dxStream << "#   mol addrep top\n";
    dxStream << "#\n";
    dxStream << "# To load a time series, use the companion .vmd script.\n";
    dxStream << "#\n";
    dxStream << "# In Jmol:\n";
    dxStream << "#   isosurface cutoff 0.85 \"" << fname << "\"\n";
    dxStream << "#\n";

    // Positions object
    dxStream << "object 1 class gridpositions counts "
             << nBinsX_ << " " << nBinsY_ << " " << nBinsZ_ << "\n";
    dxStream << "origin "
             << std::fixed << std::setprecision(6)
             << origin[0] << " " << origin[1] << " " << origin[2] << "\n";
    dxStream << "delta " << voxelSize_ << " 0 0\n";
    dxStream << "delta 0 " << voxelSize_ << " 0\n";
    dxStream << "delta 0 0 " << voxelSize_ << "\n";

    // Connections object
    dxStream << "object 2 class gridconnections counts "
             << nBinsX_ << " " << nBinsY_ << " " << nBinsZ_ << "\n";

    // Data array
    dxStream << "object 3 class array type double rank 0 items "
             << nItems << " data follows\n";

    // Data values: x slow, y medium, z fast (DX convention)
    int count = 0;
    for (int i = 0; i < nBinsX_; i++) {
      for (int j = 0; j < nBinsY_; j++) {
        for (int k = 0; k < nBinsZ_; k++) {
          dxStream << std::scientific << std::setprecision(6)
                   << field[i][j][k];
          count++;
          if (count % 3 == 0)
            dxStream << "\n";
          else
            dxStream << " ";
        }
      }
    }
    if (count % 3 != 0) dxStream << "\n";

    dxStream << "attribute \"dep\" string \"positions\"\n";
    dxStream << "object \"Tetrahedrality Q\" class field\n";
    dxStream << "component \"positions\" value 1\n";
    dxStream << "component \"connections\" value 2\n";
    dxStream << "component \"data\" value 3\n";

    dxStream.close();
  }

  void Qsurf::writeVMDScript(const std::string& scriptName,
                             const std::string& baseName,
                             int nFramesWritten) {
    std::ofstream tcl(scriptName.c_str());
    if (!tcl.is_open()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Qsurf: unable to open %s for writing\n",
               scriptName.c_str());
      painCave.isFatal = 1;
      simError();
    }

    Revision r;
    tcl << "#!/usr/bin/env vmd\n";
    tcl << "# -------------------------------------------------------\n";
    tcl << "# VMD Tcl script for animating Qsurf tetrahedrality field\n";
    tcl << "# Generated by OpenMD " << r.getFullRevision() << "\n";
    tcl << "# " << r.getBuildDate() << "\n";
    tcl << "# -------------------------------------------------------\n";
    tcl << "#\n";
    tcl << "# Usage:\n";
    tcl << "#   vmd -e " << scriptName << "\n";
    tcl << "#\n";
    tcl << "# This script loads " << nFramesWritten
        << " volumetric Q-field snapshots\n";
    tcl << "# into a single VMD molecule, creates an Isosurface\n";
    tcl << "# representation, and installs a trace callback so the\n";
    tcl << "# displayed volumetric dataset updates as you scrub\n";
    tcl << "# through frames.\n";
    tcl << "#\n";
    tcl << "# Once loaded, you can:\n";
    tcl << "#   - Adjust the isovalue in Graphics > Representations\n";
    tcl << "#   - Change Draw style (Solid/Wireframe/Points)\n";
    tcl << "#   - Color by Volume to see Q variation on the surface\n";
    tcl << "#   - Overlay the molecular trajectory by loading the\n";
    tcl << "#     .omd dump into the SAME molecule first, e.g.:\n";
    tcl << "#       mol new system.omd waitfor all\n";
    tcl << "#     then source this script with the molecule selected.\n";
    tcl << "# -------------------------------------------------------\n";
    tcl << "\n";

    // Load the first frame to create the molecule
    tcl << "# Load the first DX file as a new molecule\n";
    tcl << "set qmol [mol new {" << baseName
        << ".0000.dx} type dx waitfor all]\n\n";

    // Load remaining frames
    if (nFramesWritten > 1) {
      tcl << "# Load remaining frames into the same molecule\n";
      tcl << "for {set i 1} {$i < " << nFramesWritten << "} {incr i} {\n";
      tcl << "  set fname [format \"" << baseName
          << ".%04d.dx\" $i]\n";
      tcl << "  mol addfile $fname type dx waitfor all molid $qmol\n";
      tcl << "}\n\n";
    }

    // Create the isosurface representation
    tcl << "# Create an isosurface representation\n";
    tcl << "# Default isovalue = 0.85 (adjust to taste)\n";
    tcl << "mol delrep 0 $qmol\n";
    tcl << "mol representation Isosurface 0.85 0 0 0 1 1\n";
    tcl << "mol color Volume 0\n";
    tcl << "mol selection {all}\n";
    tcl << "mol material Opaque\n";
    tcl << "mol addrep $qmol\n\n";

    // Store the rep name for the trace callback
    tcl << "# Get the representation name for later updates\n";
    tcl << "set qrep [mol repname $qmol 0]\n\n";

    // Install the trace callback
    tcl << "# -------------------------------------------------------\n";
    tcl << "# Trace callback: when the frame changes, update the\n";
    tcl << "# isosurface to use the volumetric dataset corresponding\n";
    tcl << "# to the current frame number.\n";
    tcl << "# -------------------------------------------------------\n";
    tcl << "proc qsurf_update_frame {args} {\n";
    tcl << "  global qmol qrep\n";
    tcl << "  set f [molinfo $qmol get frame]\n";
    tcl << "  set repid [mol repindex $qmol $qrep]\n";
    tcl << "  if {$repid >= 0} {\n";
    tcl << "    mol modstyle $repid $qmol Isosurface 0.85 $f 0 0 1 1\n";
    tcl << "  }\n";
    tcl << "}\n\n";
    tcl << "trace variable vmd_frame($qmol) w qsurf_update_frame\n\n";

    tcl << "# Set animation to first frame\n";
    tcl << "animate goto 0\n";
    tcl << "qsurf_update_frame\n\n";

    tcl << "puts \"Qsurf: loaded " << nFramesWritten
        << " frames into molecule $qmol\"\n";
    tcl << "puts \"Qsurf: use the animation controls to scrub through "
        << "frames\"\n";
    tcl << "puts \"Qsurf: adjust isovalue in Graphics > Representations"
        << "\"\n";

    tcl.close();
  }

}  // namespace OpenMD
