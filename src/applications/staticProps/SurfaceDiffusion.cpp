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

/* Surface Diffusion 
 * Attempting to track/measure the surface diffusion rates of particles on... wait for it..
 * a surface.
 * This program was initially created to track Platinum particles moving around a 557 surface.
 * Hence why we are trying to keep the x and y movement separate. 
 *
 */

#include <algorithm>
#include <fstream>
#include "applications/staticProps/SurfaceDiffusion.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
namespace OpenMD {
  
  SurfaceDiffusion::SurfaceDiffusion(SimInfo* info, const std::string& filename, const std::string& sele, RealType len)
    : StaticAnalyser(info, filename, 1), selectionScript_(sele),  evaluator_(info), seleMan1_(info){

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator_.evaluate());
    }       
    
    //Depending on the selection 'sele1="select Pt"' need a vector equal to the
    //number of Platinums in the system (for this specific case)
    selectionCount_ = seleMan1_.getSelectionCount();
    cout << "SelectionCount_: " << selectionCount_ << "\n";
    
    moBool_.resize(selectionCount_);
    positions_.resize(selectionCount_);

    filename_ = filename; 
    singleMoveDistance_ = 2.0;
  }

  SurfaceDiffusion::~SurfaceDiffusion(){
  
  }

  void SurfaceDiffusion::process() {
    StuntDouble* sd;
    bool usePeriodicBoundaryConditions_ = info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    frames_ = 0;
    nProcessed_ = nFrames/step_;

    // positions_ and moBool_ are 2D arrays, need the second dimension
    // filled as well
    for(int i = 0; i < selectionCount_; i++){
      moBool_[i].resize(nFrames);
      positions_[i].resize(nFrames);
    }

    int iterator;
    int index = 0; 
    /* Loop over all frames storing the positions in a vec< vec<pos> >
     * At the end, positions.length() should equal seleMan1_.size() or
     * w/e And positions[index].length() should equal nFrames (or
     * nFrames/istep)
     */
    for(int istep = 0; istep < nFrames; istep += step_){
      frames_++;
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      
      index = 0; // count over atoms since iterators aren't the most
		 // friendly for such plebian things
      for(sd = seleMan1_.beginSelected(iterator); sd != NULL; 
	  sd = seleMan1_.nextSelected(iterator)){
	Vector3d pos = sd->getPos();
	positions_[index][istep] = pos;
	index++;
      }      
    }

    cout << "Position Array size: " << positions_.size() << "\n";
    cout << "Frames analyzed: " << positions_[0].size() << "\n";
      
    for(std::size_t i = 0; i < positions_.size(); i++){
      int frameIndex = positions_[i].size();
      for(int j = 1; j < frameIndex; j++){
	Vector3d posF1 = positions_[i][j-1];
	Vector3d posF2 = positions_[i][j];
	Vector3d diff = posF2 - posF1;
	if(usePeriodicBoundaryConditions_){
	  currentSnapshot_->wrapVector(diff);
	}
	double dist = diff.length();
	if(dist > singleMoveDistance_){
	  moBool_[i][j] = true;
	}else{
	  moBool_[i][j] = false;
	}
      }
    } 

    int mobileAtomCount = 0;
    for(std::size_t i = 0; i < moBool_.size(); i++){
      int frameIndex = moBool_[i].size();
      bool mobileAtom = false;
      for(int j = 0; j < frameIndex; j++){
	mobileAtom = mobileAtom || moBool_[i][j];
      }
      moBool_[i][0] = mobileAtom; // is true if any value later in the
				  // array is true, false otherwise
      if(mobileAtom){
	mobileAtomCount++;
      }
    } 

    cout << "Mobile atom count: " << mobileAtomCount << "\n";

    // Here I shrink the size of the arrays, why look through 3888,
    // when you only need ~800.  Additionally, all of these are mobile
    // at some point in time, the others aren't, dead weight and
    // memory
    positions2_.resize(mobileAtomCount);
    moBool2_.resize(mobileAtomCount);
    int pos2index = 0;
    for(std::size_t i = 0; i < positions_.size(); i++){
      int frameCount = positions_[i].size();
      if(moBool_[i][0]){
	for(int j = 0; j < frameCount; j++){
	  positions2_[pos2index].push_back(positions_[i][j]);
	  moBool2_[pos2index].push_back(moBool_[i][j]);
	}
	pos2index++;
      }
    }	
    
    positions_.clear();
    moBool_.clear();

    cout << "positions_ has been cleared: " << positions_.size() << "\n";
    cout << "positions2_ has been filled: " << positions2_.size() << "\n";
    cout << "positions2_ has " << positions2_[0].size() << " frames\n";
  
    //The important one! 
    positionCorrelation();
 

    //Write out my data
    std::ofstream diffStream;
    setOutputName(getPrefix(filename_) + ".Mdiffusion");
    diffStream.open(outputFilename_.c_str());
    diffStream << "#X&Y diffusion amounts\n";
    diffStream << "#singleMoveDistance_: " << singleMoveDistance_ << "\n";
    diffStream << "#Number of mobile atoms: " << positions2_.size() << "\n";
    diffStream << "#time, <x(t)-x(0)>, <y(t)-y(0)>, <r(t)-r(0)>\n";

    for(std::size_t i = 0; i < xHist_.size(); i++){
      diffStream << i << ", " << xHist_[i] << ", " << yHist_[i] << ", " 
		 << rHist_[i] << "\n";
    }
    diffStream.close();

  }
 
  void SurfaceDiffusion::positionCorrelation(){
    RealType xDist = 0.0;
    RealType yDist = 0.0;
    RealType rDist = 0.0;
    int timeShift = 0;
    Vector3d kPos;
    Vector3d jPos;
    //biggest timeShift is positions2_[0].size() - 1?
    xHist_.clear();
    yHist_.clear();
    rHist_.clear();
    count_.clear();
    int frameResize = positions2_[0].size();
    xHist_.resize(frameResize);
    yHist_.resize(frameResize);
    rHist_.resize(frameResize);
    count_.resize(frameResize);
    //loop over particles
    // loop over frames starting at j
    //  loop over frames starting at k = j (time shift of 0)
    for(std::size_t i = 0; i < positions2_.size(); i++){
      int frames = positions2_[i].size() - 1; // for counting
					      // properly, otherwise
					      // moBool2_[i][j+1] will
					      // go over
      for(int j = 0; j < frames; j++){
	// if the particle is mobile between j and j + 1, then count
	// it for all timeShifts
	if(moBool2_[i][j+1]){
	  for(std::size_t k = j; k < positions2_[0].size(); k++){
	    //<x(t)-x(0)>  <y(t)-y(0)>  <r(t)-r(0)>
	    //The positions stored are not wrapped, thus I don't need
	    //to worry about pbc	  
	    //Mean square displacement
	    //So I do want the squared distances
	  
	    kPos = positions2_[i][k];
	    jPos = positions2_[i][j];
	    xDist = kPos.x() - jPos.x();
	    xDist = xDist*xDist;
	  
	    yDist = kPos.y() - jPos.y();
	    yDist = yDist*yDist;

	    rDist = (kPos - jPos).lengthSquare();
  

	    timeShift = k - j;
	    xHist_[timeShift] += xDist;
	    yHist_[timeShift] += yDist;  
	    rHist_[timeShift] += rDist;
	    count_[timeShift]++;
	  }
	}
      }
    }
    cout << "X, Y, R calculated\n";

    for(std::size_t i = 0; i < xHist_.size(); i++){
      xHist_[i] = xHist_[i]/(count_[i]);
      yHist_[i] = yHist_[i]/(count_[i]);
      rHist_[i] = rHist_[i]/(count_[i]);
    }
    cout << "X, Y, R normalized\n";
  }

}
