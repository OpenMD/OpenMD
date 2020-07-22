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
 
#include "applications/staticProps/TetrahedralityHBMatrix.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Revision.hpp"
#include "utils/Constants.hpp"
#include <vector>

namespace OpenMD {

  TetrahedralityHBMatrix::TetrahedralityHBMatrix(SimInfo* info, 
                                                 const std::string& filename, 
                                                 const std::string& sele,
                                                 double rCut, double OOcut,
                                                 double thetaCut, double OHcut,
                                                 int nbins) : 
    StaticAnalyser(info, filename, nbins), selectionScript_(sele), 
    seleMan_(info), evaluator_(info), rCut_(rCut),  OOCut_(OOcut),
    thetaCut_(thetaCut), OHCut_(OHcut) {

    setAnalysisType("Tetrahedrality HBond Matrix");   
    setOutputName(getPrefix(filename) + ".hbq");
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    Q_histogram_.resize(nBins_);

    for (unsigned int i = 0; i < nBins_; i++) {
      Q_histogram_[i].resize(nBins_);
    }

    // Q can take values from 0 to 1

    MinQ_ = 0.0;
    MaxQ_ = 1.1;
    deltaQ_ = (MaxQ_ - MinQ_) / nBins_;

    std::stringstream params;
    params << " rCut = " << rCut_
           << ", OOcut = " << OOCut_
           << ", thetacut = " << thetaCut_
           << ", OHcut = " << OHCut_
           << ", nbins = " << nBins_
           << ", deltaQ = " << deltaQ_;
    const std::string paramString = params.str();
    setParameterString( paramString );        
  }

  TetrahedralityHBMatrix::~TetrahedralityHBMatrix() {
    for (unsigned int i = 0; i < nBins_; i++) {
      Q_histogram_[i].clear();
    }    
    Q_histogram_.clear(); 
  }
  
  void TetrahedralityHBMatrix::initializeHistogram() {
    for (unsigned int i = 0; i < nBins_; i++) {
      std::fill(Q_histogram_[i].begin(), Q_histogram_[i].end(), 0);
    }
  }
  
  void TetrahedralityHBMatrix::process() {
    Molecule* mol1;
    Molecule* mol2;
    Molecule* moli;
    Molecule* molj;
    SimInfo::MoleculeIterator mi;
    Vector3d vec;
    Vector3d ri, rj, rk, rik, rkj, dposition, tposition;
    RealType r, cospsi;
    std::vector<Molecule::HBondDonor*>::iterator hbdi;
    Molecule::HBondDonor* hbd;
    std::vector<Atom*>::iterator hbai;
    Atom* hba;
    Vector3d dPos;
    Vector3d aPos;
    Vector3d hPos;
    Vector3d DH;
    Vector3d DA;
    Vector3d HA;
    Vector3d uDA;
    RealType DAdist, DHdist, HAdist, theta, ctheta, Qk, q1, q2;
    std::vector<std::pair<RealType,Molecule*> > myNeighbors;
    int myIndex, ii, jj, index1, index2;
    
    bool usePeriodicBoundaryConditions_ = info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      Q_.clear();
      Q_.resize( info_->getNGlobalMolecules(), 0.0 );   
      
      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      // outer loop is over the selected Molecules:
      
      for (mol1 = seleMan_.beginSelectedMolecule(ii);
           mol1 != NULL; mol1 = seleMan_.nextSelectedMolecule(ii)) {
	
        myIndex = mol1->getGlobalIndex();
        
	Qk = 1.0;
        
	myNeighbors.clear();
        
        // inner loop is over all Molecules in the system:
        
        for (mol2 = info_->beginMolecule(mi); mol2 != NULL; 
             mol2 = info_->nextMolecule(mi)) {
          
          if (mol2->getGlobalIndex() != myIndex) {
            
            vec = mol1->getCom() - mol2->getCom();       
	    
            if (usePeriodicBoundaryConditions_) 
              currentSnapshot_->wrapVector(vec);
            
            r = vec.length();             
            
            // Check to see if neighbor is in bond cutoff 
            
            if (r < rCut_) { 
              
              myNeighbors.push_back(std::make_pair(r,mol2));
            }
          }
        }      

	// Sort the vector using predicate and std::sort
	std::sort(myNeighbors.begin(), myNeighbors.end());
	
	// Use only the 4 closest neighbors to do the rest of the work:	
	// int nbors =  myNeighbors.size()> 4 ? 4 : myNeighbors.size();
	int nbors = myNeighbors.size();
	int nang = int (0.5 * (nbors * (nbors - 1)));

	rk = mol1->getCom();

	for (int i = 0; i < nbors-1; i++) {	  

	  moli = myNeighbors[i].second;
	  ri = moli->getCom();
	  rik = rk - ri;
	  if (usePeriodicBoundaryConditions_) 
	    currentSnapshot_->wrapVector(rik);
	  
	  rik.normalize();

	  for (int j = i+1; j < nbors; j++) {	    

	    molj = myNeighbors[j].second;
	    rj = molj->getCom();
	    rkj = rk - rj;
	    if (usePeriodicBoundaryConditions_) 
	      currentSnapshot_->wrapVector(rkj);
	    rkj.normalize();
	    
	    cospsi = dot(rik,rkj);
            
	    // Calculates scaled Qk for each molecule using calculated
	    // angles from the actual number of nearest neighbors.
	    Qk = Qk - (pow(cospsi + 1.0 / 3.0, 2) * 9.0 / (4.0 * nang) );
	  }
	}
        
	if (nang > 0) {
          Q_[myIndex] = Qk;
          
        }
      }
      // Now to scan for histogram:

      for (mol1 = seleMan_.beginSelectedMolecule(ii);
           mol1 != NULL; mol1 = seleMan_.nextSelectedMolecule(ii)) {
        
        for (mol2 = seleMan_.beginSelectedMolecule(jj);
             mol2 != NULL; mol2 = seleMan_.nextSelectedMolecule(jj)) {
          
          // loop over the possible donors in molecule 1:
          for (hbd = mol1->beginHBondDonor(hbdi); hbd != NULL;
               hbd = mol1->nextHBondDonor(hbdi)) {
            dPos = hbd->donorAtom->getPos(); 
            hPos = hbd->donatedHydrogen->getPos();
            DH = hPos - dPos; 
            currentSnapshot_->wrapVector(DH);
            DHdist = DH.length();
            
            // loop over the possible acceptors in molecule 2:
            for (hba = mol2->beginHBondAcceptor(hbai); hba != NULL;
                 hba = mol2->nextHBondAcceptor(hbai)) {
              aPos = hba->getPos();
              DA = aPos - dPos;
              currentSnapshot_->wrapVector(DA);
              DAdist = DA.length();
              
              // Distance criteria: are the donor and acceptor atoms
              // close enough?
              if (DAdist < OOCut_) {
                HA = aPos - hPos;
                currentSnapshot_->wrapVector(HA);
                HAdist = HA.length();
                
                ctheta = dot(DH, DA) / (DHdist * DAdist);
                theta = acos(ctheta) * 180.0 / Constants::PI;
              
                // Angle criteria: are the D-H and D-A and vectors close?
                if (theta < thetaCut_ && HAdist < OHCut_) {
                  // We have a hydrogen bond!                
                  index1 = mol1->getGlobalIndex();
                  index2 = mol2->getGlobalIndex();
                  q1 = Q_[index1];
                  q2 = Q_[index2];
                  collectHistogram(q1, q2);                
                }
              }
            }
          }
          // now loop over the possible acceptors in molecule 1:
          for (hba = mol1->beginHBondAcceptor(hbai); hba != NULL;
               hba = mol1->nextHBondAcceptor(hbai)) {
            aPos = hba->getPos();
            
            // loop over the possible donors in molecule 2:
            for (hbd = mol2->beginHBondDonor(hbdi); hbd != NULL;
                 hbd = mol2->nextHBondDonor(hbdi)) {
              dPos = hbd->donorAtom->getPos();
              
              DA = aPos - dPos;
              currentSnapshot_->wrapVector(DA);
              DAdist = DA.length();
              
              // Distance criteria: are the donor and acceptor atoms
              // close enough?
              if (DAdist < OOCut_) {
                hPos = hbd->donatedHydrogen->getPos();
                HA = aPos - hPos;
                currentSnapshot_->wrapVector(HA);
                HAdist = HA.length();
                
                DH = hPos - dPos; 
                currentSnapshot_->wrapVector(DH);
                DHdist = DH.length();
                ctheta = dot(DH, DA) / (DHdist * DAdist);
                theta = acos(ctheta) * 180.0 / Constants::PI;              
                // Angle criteria: are the D-H and D-A and vectors close?
                if (theta < thetaCut_ && HAdist < OHCut_) {
                  // We have a hydrogen bond!
                  // keep these indexed by donor then acceptor:
                  index1 = mol2->getGlobalIndex();
                  index2 = mol1->getGlobalIndex();
                  q1 = Q_[index1];
                  q2 = Q_[index2];
                  collectHistogram(q1, q2);
                }
              }
            }
          }
        }           
      }            
    }
    writeOutput();
  }
        
  void TetrahedralityHBMatrix::collectHistogram(RealType q1, RealType q2) {

    if (q1 > MinQ_ && q1 < MaxQ_) {     
      int bin1 = int((q1 - MinQ_) / deltaQ_);
      if (q2 > MinQ_ && q2 < MaxQ_) {     
        int bin2 = int((q2 - MinQ_) / deltaQ_);
        Q_histogram_[bin1][bin2] += 1;
        count_++;
      } else {
        // cerr << "q2 = " << q2 << "\n";
      }
    } else {
      // cerr << "q1 = " << q1 << "\n";
    }
  }

  void TetrahedralityHBMatrix::writeOutput() {
    std::ofstream ofs(outputFilename_.c_str());
    if (ofs.is_open()) {
      Revision r;
      ofs << "# " << getAnalysisType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script: \"" << selectionScript_ << "\"\n";
      
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";

      ofs << "# Hydrogen Bond count: " << count_ << "\n";
      
      for (unsigned int i = 0; i < Q_histogram_.size(); ++i) {
        for(unsigned int j = 0; j < Q_histogram_[i].size(); ++j) {
          ofs <<  RealType(Q_histogram_[i][j])/RealType(count_) << "\t";
        }
        ofs << "\n";        
      }
      
    } else {
      sprintf(painCave.errMsg, "TetrahedralityHBMatrix: unable to open %s\n", 
              outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }
    
    ofs.close();
  }
}



