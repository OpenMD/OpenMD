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
 
#include "applications/staticProps/TetrahedralityParamXYZ.hpp"
#include "utils/simError.h"
#include "utils/Constants.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;
namespace OpenMD {
  TetrahedralityParamXYZ::TetrahedralityParamXYZ(SimInfo* info,  
                                                 const std::string& filename, 
                                                 const std::string& sele1,
                                                 const std::string& sele2,
                                                 RealType rCut, 
                                                 RealType voxelSize,
                                                 RealType gaussWidth) 
    : StaticAnalyser(info, filename, 1), 
      selectionScript1_(sele1), selectionScript2_(sele2), 
      seleMan1_(info),  seleMan2_(info), evaluator1_(info), evaluator2_(info),
      rCut_(rCut), voxelSize_(voxelSize), gaussWidth_(gaussWidth) {
    
    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }
    evaluator2_.loadScriptString(sele2);
    if (!evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    }
    
    Mat3x3d hmat = info->getSnapshotManager()->getCurrentSnapshot()->getHmat();

    nBins_(0) = int(hmat(0,0) / voxelSize);
    nBins_(1) = int(hmat(1,1) / voxelSize);
    nBins_(2) = int(hmat(2,2) / voxelSize);

    hist_.resize(nBins_(0));
    count_.resize(nBins_(0));
    for (int i = 0 ; i < nBins_(0); ++i) {
      hist_[i].resize(nBins_(1));
      count_[i].resize(nBins_(1));
      for(int j = 0; j < nBins_(1); ++j) {
        hist_[i][j].resize(nBins_(2));
        count_[i][j].resize(nBins_(2));
        std::fill(hist_[i][j].begin(), hist_[i][j].end(), 0.0);
        std::fill(count_[i][j].begin(), count_[i][j].end(), 0.0);

      }
    }   
    
    setOutputName(getPrefix(filename) + ".Qxyz");
  }
  
  TetrahedralityParamXYZ::~TetrahedralityParamXYZ() {
  }
    
  void TetrahedralityParamXYZ::process() {
    StuntDouble* sd;
    StuntDouble* sd2;
    StuntDouble* sdi;
    StuntDouble* sdj;
    int myIndex;
    Vector3d vec;
    Vector3d ri, rj, rk, rik, rkj;
    RealType r;
    RealType cospsi;
    RealType Qk;
    std::vector<std::pair<RealType,StuntDouble*> > myNeighbors;
    //std::vector<std::pair<Vector3d, RealType> > qvals;
    //std::vector<std::pair<Vector3d, RealType> >::iterator qiter;
    int isd1;
    int isd2;
    bool usePeriodicBoundaryConditions_ = info_->getSimParams()->getUsePeriodicBoundaryConditions();


    int kMax = int(5.0 * gaussWidth_ / voxelSize_);
    int kSqLim = kMax*kMax;
    cerr << "gw = " << gaussWidth_ << " vS = " << voxelSize_ << " kMax = " 
	 << kMax << " kSqLim = " << kSqLim << "\n";
    
    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);

      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      Mat3x3d hmat = currentSnapshot_->getHmat();
      Vector3d halfBox = Vector3d(hmat(0,0), hmat(1,1), hmat(2,2)) / 2.0;

      if (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
      }
      
      if (evaluator2_.isDynamic()) {
        seleMan2_.setSelectionSet(evaluator2_.evaluate());
      }
      
      //qvals.clear();

      // outer loop is over the selected StuntDoubles:
      for (sd = seleMan1_.beginSelected(isd1); sd != NULL;
           sd = seleMan1_.nextSelected(isd1)) {
        
        myIndex = sd->getGlobalIndex();
        
        Qk = 1.0;	  
        myNeighbors.clear();       

        for (sd2 = seleMan2_.beginSelected(isd2); sd2 != NULL;
             sd2 = seleMan2_.nextSelected(isd2)) {
          
          if (sd2->getGlobalIndex() != myIndex) {
            
            vec = sd->getPos() - sd2->getPos();       
            
            if (usePeriodicBoundaryConditions_) 
              currentSnapshot_->wrapVector(vec);
            
            r = vec.length();             
            
            // Check to see if neighbor is in bond cutoff 
            
            if (r < rCut_) {                
              myNeighbors.push_back(std::make_pair(r,sd2));
            }
          }
        }
        
        // Sort the vector using predicate and std::sort
        std::sort(myNeighbors.begin(), myNeighbors.end());
        
        // Use only the 4 closest neighbors to do the rest of the work:
        
        int nbors =  myNeighbors.size()> 4 ? 4 : myNeighbors.size();
        int nang = int (0.5 * (nbors * (nbors - 1)));
        
        rk = sd->getPos();
        
        for (int i = 0; i < nbors-1; i++) {       
          
          sdi = myNeighbors[i].second;
          ri = sdi->getPos();
          rik = rk - ri;
          if (usePeriodicBoundaryConditions_) 
            currentSnapshot_->wrapVector(rik);
          
          rik.normalize();
          
          for (int j = i+1; j < nbors; j++) {       
            
            sdj = myNeighbors[j].second;
            rj = sdj->getPos();
            rkj = rk - rj;
            if (usePeriodicBoundaryConditions_) 
              currentSnapshot_->wrapVector(rkj);
            rkj.normalize();
            
            cospsi = dot(rik,rkj);           
            
            // Calculates scaled Qk for each molecule using calculated
            // angles from 4 or fewer nearest neighbors.
            Qk -=  (pow(cospsi + 1.0 / 3.0, 2) * 2.25 / nang);            
          }
        }
        
        if (nang > 0) {
          if (usePeriodicBoundaryConditions_)
            currentSnapshot_->wrapVector(rk);
          //qvals.push_back(std::make_pair(rk, Qk));
          
          Vector3d pos = rk + halfBox;


          Vector3i whichVoxel(int(pos[0] / voxelSize_), 
                              int(pos[1] / voxelSize_), 
                              int(pos[2] / voxelSize_));
          
          for (int l = -kMax; l <= kMax; l++) {
            for (int m = -kMax; m <= kMax; m++) {
              for (int n = -kMax; n <= kMax; n++) {
                int kk = l*l + m*m + n*n;
                if(kk <= kSqLim) {

                  int ll = (whichVoxel[0] + l) % nBins_(0);
                  ll = ll < 0 ? nBins_(0) + ll : ll;
                  int mm = (whichVoxel[1] + m) % nBins_(1);
                  mm = mm < 0 ? nBins_(1) + mm : mm;
                  int nn = (whichVoxel[2] + n) % nBins_(2);
                  nn = nn < 0 ? nBins_(2) + nn : nn;

                  Vector3d bPos = Vector3d(ll,mm,nn) * voxelSize_ - halfBox;
                  Vector3d d = bPos - rk;
                  currentSnapshot_->wrapVector(d);
                  RealType denom = pow(2.0 * sqrt(Constants::PI) * gaussWidth_, 3);
                  RealType exponent = -dot(d,d) / pow(2.0*gaussWidth_, 2);
                  RealType weight = exp(exponent) / denom;
                  count_[ll][mm][nn] += weight;
                  hist_[ll][mm][nn] += weight * Qk;
                }
              }
            }
          }
        }
      }

      // for (int i = 0; i < nBins_(0); ++i) {
      //   for(int j = 0; j < nBins_(1); ++j) {
      //     for(int k = 0; k < nBins_(2); ++k) {
      //       Vector3d pos = Vector3d(i, j, k) * voxelSize_ - halfBox;
      //       for(qiter = qvals.begin(); qiter != qvals.end(); ++qiter) {
      //         Vector3d d = pos - (*qiter).first;
      //         currentSnapshot_->wrapVector(d);
      //         RealType denom = pow(2.0 * sqrt(Constants::PI) * gaussWidth_, 3);
      //         RealType exponent = -dot(d,d) / pow(2.0*gaussWidth_, 2);
      //         RealType weight = exp(exponent) / denom;
      //         count_[i][j][k] += weight;
      //         hist_[i][j][k] += weight * (*qiter).second;
      //       }
      //     }
      //   }
      // }
    }
    writeQxyz();
  }
  
  void TetrahedralityParamXYZ::writeQxyz() {

    Mat3x3d hmat = info_->getSnapshotManager()->getCurrentSnapshot()->getHmat();

    // normalize by total weight in voxel:
    for (unsigned int i = 0; i < hist_.size(); ++i) { 
      for(unsigned int j = 0; j < hist_[i].size(); ++j) { 
        for(unsigned int k = 0;k < hist_[i][j].size(); ++k) {
          hist_[i][j][k] = hist_[i][j][k] / count_[i][j][k];
        }
      }
    }

    std::ofstream qXYZstream(outputFilename_.c_str());
    if (qXYZstream.is_open()) {
      qXYZstream <<  "# AmiraMesh ASCII 1.0\n\n";
      qXYZstream <<  "# Dimensions in x-, y-, and z-direction\n";
      qXYZstream <<  "  define Lattice " << hist_.size() << " " 
		 << hist_[0].size() << " " << hist_[0][0].size() << "\n";
      
      qXYZstream <<  "Parameters {\n";
      qXYZstream <<  "    CoordType \"uniform\",\n";
      qXYZstream <<  "    # BoundingBox is xmin xmax ymin ymax zmin zmax\n";
      qXYZstream <<  "    BoundingBox 0.0 " << hmat(0,0) << 
	" 0.0 " << hmat(1,1) << 
	" 0.0 " << hmat(2,2) << "\n";
      qXYZstream <<  "}\n";
      
      qXYZstream <<  "Lattice { double ScalarField } = @1\n";
      
      qXYZstream <<  "@1\n";

      for (std::size_t k = 0; k < hist_[0][0].size(); ++k) { 
        for(std::size_t j = 0; j < hist_[0].size(); ++j) { 
          for(std::size_t i = 0; i < hist_.size(); ++i) {
	    qXYZstream << hist_[i][j][k] << " ";

            //qXYZstream.write(reinterpret_cast<char *>( &hist_[i][j][k] ),
	    // sizeof( hist_[i][j][k] ));
          }
        }
      }
      
    } else {      
      sprintf(painCave.errMsg, "TetrahedralityParamXYZ: unable to open %s\n", 
              outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }    
    qXYZstream.close();
  }
}



