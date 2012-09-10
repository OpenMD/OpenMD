/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [4] , Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). *
 *  Created by J. Daniel Gezelter on 09/26/06.
 *  @author  J. Daniel Gezelter
 *  @version $Id: BondOrderParameter.cpp 1442 2010-05-10 17:28:26Z gezelter $
 *
 */
 
#include "applications/staticProps/TetrahedralityParamZ.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/NumericConstant.hpp"
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;

namespace OpenMD 
{
  TetrahedralityParamZ::TetrahedralityParamZ(SimInfo* info, 
                                           const std::string& filename, 
                                           const std::string& sele,
                                           double rCut, int nzbins) : StaticAnalyser(info, filename), selectionScript_(sele), evaluator_(info), seleMan1_(info),seleMan2_(info), nZBins_(nzbins)
  {
    //nZBins_ = 50;
    //std ::cerrnZBins_:"<<nZBins_<<"\t"<<"nzbins:"<<nzbins<<endl;
    // nZBins_ = 90;
    //fixed numbe of bins
    count_.resize(nZBins_);
    sliceSDLists_.resize(nZBins_);
    Qave_.resize(nZBins_);

    setOutputName(getPrefix(filename) + ".q");
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) 
      {
	seleMan1_.setSelectionSet(evaluator_.evaluate());
	seleMan2_.setSelectionSet(evaluator_.evaluate());
      }

    // Set up cutoff radius:
    rCut_ = rCut;

    // Q can take values from 0 to 1
    MinQ_ = 0.0;
    MaxQ_ = 1.1;
    deltaQ_ = (MaxQ_ - MinQ_)/nzbins;
  }

  TetrahedralityParamZ::~TetrahedralityParamZ() 
  {
    Q_histogram_.clear(); 
  }
  
  void TetrahedralityParamZ::initalizeHistogram() 
  {
    std::fill(Q_histogram_.begin(), Q_histogram_.end(), 0);
  }
  

  

  void TetrahedralityParamZ::process() 
  {
    Molecule* mol;
    StuntDouble* sd;
    StuntDouble* sd2;
    StuntDouble* sdi;
    StuntDouble* sdj;
    RigidBody* rb;
    int myIndex;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    Vector3d vec;
    Vector3d ri, rj, rk, rik, rkj, dposition, tposition;
    RealType r;
    RealType cospsi;
    RealType Qk;

    std::vector<std::pair<RealType,StuntDouble*> > myNeighbors;
    int isd1, isd2;
    cerr << "After Creation of variables in TP:process()\n";
    DumpReader reader(info_, dumpFilename_);    
    cerr << "The DumpReader was created?\n";
    cerr << "nZbins: " << nZBins_ << "\n";
    int nFrames = reader.getNFrames();
    frameCounter_ = 0;
    nProcessed_=nFrames/step_;
    reader.readFrame(0); 
    currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
    Mat3x3d hmat = currentSnapshot_->getHmat();
    zBox_.push_back(hmat(2,2));
    
    RealType halfBoxZ_ = hmat(2,2) / 2.0; 

    Distorted_.clear();
    Tetrahedral_.clear();
    int i;
    for(i=0;i<nZBins_;i++) {
      sliceSDLists_[i].clear();
    }

    //LOOP OVER ALL FRAMES 
    for (int istep = 0; istep < nFrames; istep += step_) {
      int i;
      for(i=0;i<nZBins_;i++) {
        count_[i]=0;
      }	
	
      reader.readFrame(istep);
      frameCounter_++;
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      
      if (evaluator_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator_.evaluate());
        seleMan2_.setSelectionSet(evaluator_.evaluate());
      }
	
      // update the positions of atoms which belong to the rigidbodies
      for (mol = info_->beginMolecule(mi); mol != NULL;
           mol = info_->nextMolecule(mi)) {
        for (rb = mol->beginRigidBody(rbIter); rb != NULL;
             rb = mol->nextRigidBody(rbIter)) {
          rb->updateAtoms();
        }        
      }           

      // outer loop is over the selected StuntDoubles:
      int idk=0;
      for (sd = seleMan1_.beginSelected(isd1); sd != NULL;
           sd = seleMan1_.nextSelected(isd1)) {
        myIndex = sd->getGlobalIndex();
        Qk = 1.0;	  
        myNeighbors.clear();
        for(sd2 = seleMan2_.beginSelected(isd2); sd2 != NULL; 
            sd2 = seleMan2_.nextSelected(isd2)){
          if(sd2->getGlobalIndex() != myIndex){
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
        //std::cerr << myNeighbors.size() <<  " neighbors within " << rCut_  << " A" << " \n";
        // Use only the 4 closest neighbors to do the rest of the work:	  
        int nbors =  myNeighbors.size() > 4 ? 4 : myNeighbors.size();
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
            Qk = Qk - (pow(cospsi + 1.0 / 3.0, 2) * 2.25 / nang);
          }
        }
	
        //std::cerr<<nbors<<endl;
        if (nang > 0) {
          //collectHistogram(Qk);
          
          // Saves positions of StuntDoubles & neighbors with
          // distorted coordination (low Qk value)
          if ((Qk < 0.55) && (Qk > 0.45)) {
            Distorted_.push_back(sd);
            dposition = sd->getPos();
          }
          
          // Saves positions of StuntDoubles & neighbors with
          // tetrahedral coordination (high Qk value)
          if (Qk > 0) {
            Tetrahedral_.push_back(sd);
            tposition = sd->getPos();
          }
          
        }
	
        //wrap the stuntdoubles into a cell      
        Vector3d pos = sd->getPos();
        if (usePeriodicBoundaryConditions_)
          currentSnapshot_->wrapVector(pos);
        sd->setPos(pos);
        // shift molecules by half a box to have bins start at 0
        int binNo = int(nZBins_ * (halfBoxZ_ + pos.z()) / hmat(2,2));
        // Patrick took out the "halfBoxZ_" part in the line above to below
        // int binNo = int(nZBins_ * (pos.z()) / hmat(2,2));
        sliceSDLists_[binNo].push_back(Qk);
        idk++;
      }//outer sd loop
    }//istep loop
 
    //Averaging the value of Qk in each bin
    for(int i=0; i< nZBins_; i++) {
      RealType Qsum=0;
      for (unsigned int k = 0; k < sliceSDLists_[i].size(); ++k) {	  
        Qsum=Qsum+sliceSDLists_[i][k];
        count_[i]++;
      }
      //std::cerr<<"past averagin Qk"<<endl;
      //std::cerr<<Qsum<<endl;
      if(count_[i]!=0) {
        Qave_.push_back(Qsum/count_[i]);
      }
      //std::cerr<<count[i]<<endl;
    }
    //std::cerr<<"nZBins_ = "<< nZBins_<<endl;
    //Writing bin#:<Qk> to a file
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      //rdfStream << "#QkZ\n";
      //rdfStream << "#nFrames:\t" << nProcessed_ << "\n";
      //rdfStream << "#selection: (" << selectionScript_ << ")\n";
      //rdfStream << "#z\tdensity\n";
      for (int i = 0; i < nZBins_; ++i) {
        if(count_[i]!=0) {
          rdfStream << ((hmat(2,2)*i)/nZBins_)+(hmat(2,2)/(2*nZBins_)) 
                    << "\t" << Qave_[i] << "\n";
        }
      }
    }
    
    writeOrderParameter();
    std::cerr << "number of distorted StuntDoubles = " 
              << Distorted_.size() << "\n";
    std::cerr << "number of tetrahedral StuntDoubles = " 
              << Tetrahedral_.size() << "\n";
    collectHistogram(Qk);

  }//void TetrahedralityParam::process() loop
  
  void TetrahedralityParamZ::collectHistogram(RealType Qk) {
    //if (Qk > MinQ_ && Qk < MaxQ_) 
    //  {
    //	int whichBin = int((Qk - MinQ_) / deltaQ_);
    //	Q_histogram_[whichBin] += 1;
    //  }
  }    

  void TetrahedralityParamZ::writeOrderParameter() {
    int nSelected = 0;
    std::ofstream osq((getOutputFileName() + "Q").c_str());
    if (osq.is_open()) {
      osq << "# Tetrahedrality Parameters\n";
      osq << "# selection: (" << selectionScript_ << ")\n";
      osq << "# \n";
      osq.close();
    } else {
      sprintf(painCave.errMsg, "TetrahedralityParamZ: unable to open %s\n", 
              (getOutputFileName() + "q").c_str());
      painCave.isFatal = 1;
      simError();  
    }
    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    if (nFrames == 1) {
      std::vector<StuntDouble*>::iterator iter;
      std::ofstream osd((getOutputFileName() + "dxyz").c_str());
      if (osd.is_open()) {
        osd << Distorted_.size() << "\n\n";
        
        for (iter = Distorted_.begin(); iter != Distorted_.end(); ++iter) {
          Vector3d position;
          position = (*iter)->getPos();
          osd << "O  " << "\t";
          for (unsigned int z=0; z<position.size(); z++) {
            osd << position[z] << "  " << "\t";
          }
          osd << "\n";
        }
        osd.close();
      }
      std::ofstream ost((getOutputFileName() + "txyz").c_str());
      if (ost.is_open()) {
        ost << Tetrahedral_.size() << "\n\n";	    
        for (iter = Tetrahedral_.begin(); iter != Tetrahedral_.end(); ++iter) {
          Vector3d position;		
          position = (*iter)->getPos();
          ost << "O  " << "\t";
          for (unsigned int z=0; z<position.size(); z++) {
            ost << position[z] << "  " << "\t";
          } 
          ost << "\n";
        }
        ost.close();
      }
    }
  }
}



