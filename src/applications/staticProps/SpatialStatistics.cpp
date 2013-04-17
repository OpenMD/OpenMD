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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). 
 */

#include "applications/staticProps/SpatialStatistics.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"

namespace OpenMD {
  
  SpatialStatistics::SpatialStatistics(SimInfo* info, const string& filename, 
                                       const string& sele, int nbins)
    : StaticAnalyser(info, filename), selectionScript_(sele),  evaluator_(info),
      seleMan_(info), nBins_(nbins){
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
    
    // Pre-load an OutputData for the count of objects:
    counts_ = new OutputData;
    counts_->units =  "objects";
    counts_->title =  "Objects";
    counts_->dataType = odtReal;
    counts_->dataHandling = odhTotal;
    counts_->accumulator.reserve(nBins_);
    for (int i = 0; i < nBins_; i++) 
      counts_->accumulator.push_back( new Accumulator() );
    
    setOutputName(getPrefix(filename) + ".spst");
  }

  SpatialStatistics::~SpatialStatistics() {
    vector<OutputData*>::iterator i;
    OutputData* outputData;
    
    for(outputData = beginOutputData(i); outputData; 
        outputData = nextOutputData(i)) {
      delete outputData;
    }
    data_.clear();

    delete counts_;
  }


  void SpatialStatistics::process() {

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames/step_;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      processFrame(istep);
    }
    writeOutput();
  }


  void SpatialStatistics::processFrame(int istep) {
    Molecule* mol;
    RigidBody* rb;
    StuntDouble* sd;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    int i;
    
    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) {
      
      // change the positions of atoms which belong to the rigidbodies
      
      for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
           rb = mol->nextRigidBody(rbIter)) {
        rb->updateAtoms();
      }
    }
    
    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
    
    // loop over the selected atoms:
    
    for (sd = seleMan_.beginSelected(i); sd != NULL; 
         sd = seleMan_.nextSelected(i)) {
      
      // figure out where that object is:
      
      Vector3d pos = sd->getPos();
      
      int bin = getBin(pos);
      
      // forward the work of statistics on to the subclass:
      
      processStuntDouble( sd, bin );

      dynamic_cast<Accumulator *>(counts_->accumulator[bin])->add(1);
    }      
  }
  

  void SpatialStatistics::writeOutput() {
    
    vector<OutputData*>::iterator i;
    OutputData* outputData;
    
    ofstream outStream(outputFilename_.c_str());
    if (outStream.is_open()) {
      
      //write title
      outStream << "# SPATIAL STATISTICS\n";
      outStream << "#";
      
      for(outputData = beginOutputData(i); outputData; 
          outputData = nextOutputData(i)) {
        outStream << "\t" << outputData->title << 
          "(" << outputData->units << ")";
        // add some extra tabs for column alignment
        if (outputData->dataType == odtVector3) outStream << "\t\t";
      }
      
      outStream << std::endl;
      
      outStream.precision(8);
      
      for (int j = 0; j < nBins_; j++) {        
        
        int counts = counts_->accumulator[j]->count();

        if (counts > 0) {
          for(outputData = beginOutputData(i); outputData; 
              outputData = nextOutputData(i)) {
            
            int n = outputData->accumulator[j]->count();
            if (n != 0) {
              writeData( outStream, outputData, j );
            }
          }
          outStream << std::endl;
        }
      }
        
      outStream << "#######################################################\n";
      outStream << "# Standard Deviations in those quantities follow:\n";
      outStream << "#######################################################\n";
      
      for (int j = 0; j < nBins_; j++) {
        int counts = counts_->accumulator[j]->count();
        if (counts > 0) {
          
          outStream << "#";
          for(outputData = beginOutputData(i); outputData; 
              outputData = nextOutputData(i)) {
            
            int n = outputData->accumulator[j]->count();
            if (n != 0) {
              writeStdDev( outStream, outputData, j );
            }
          }
          outStream << std::endl;
        }
      }
      
      outStream.flush();
      outStream.close();      
      
    } else {      
      sprintf(painCave.errMsg, "SpatialStatistics: unable to open %s\n", 
              outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }   
  }
  
  
  void SpatialStatistics::writeData(ostream& os, OutputData* dat, 
                                    unsigned int bin) {
    assert(int(bin) < nBins_);
    int n = dat->accumulator[bin]->count();
    if (n == 0) return;

    RealType r;
    Vector3d v;

    if( dat->dataType == odtReal ) {
      dynamic_cast<Accumulator*>(dat->accumulator[bin])->getAverage(r);      
      if (isinf(r) || isnan(r) ) {      
        sprintf( painCave.errMsg,
                 "SpatialStatistics detected a numerical error writing:\n"
                 "\t%s for bin %d",
                 dat->title.c_str(), bin);
        painCave.isFatal = 1;
        simError();
      }
      if (dat->dataHandling == odhTotal) r *= dat->accumulator[bin]->count();
      os << "\t" << r;      

    } else if ( dat->dataType == odtVector3 ) {
      dynamic_cast<VectorAccumulator*>(dat->accumulator[bin])->getAverage(v);
      if (isinf(v[0]) || isnan(v[0]) || 
          isinf(v[1]) || isnan(v[1]) || 
          isinf(v[2]) || isnan(v[2]) ) {      
        sprintf( painCave.errMsg,
                 "SpatialStatistics detected a numerical error writing:\n"
                 "\t%s for bin %d",
                 dat->title.c_str(), bin);
        painCave.isFatal = 1;
        simError();
      }
      if (dat->dataHandling == odhTotal) v *= dat->accumulator[bin]->count();
      os << "\t" << v[0] << "\t" << v[1] << "\t" << v[2];
    }
  }

  void SpatialStatistics::writeStdDev(ostream& os, OutputData* dat, 
                                    unsigned int bin) {
    assert(int(bin) < nBins_);
    int n = dat->accumulator[bin]->count();
    if (n == 0) return;

    RealType r;
    Vector3d v;

    if( dat->dataType == odtReal ) {
      dynamic_cast<Accumulator*>(dat->accumulator[bin])->getStdDev(r);      
      if (isinf(r) || isnan(r) ) {      
        sprintf( painCave.errMsg,
                 "SpatialStatistics detected a numerical error writing:\n"
                 "\tstandard deviation of %s for bin %d",
                 dat->title.c_str(), bin);
        painCave.isFatal = 1;
        simError();
      }
      if (dat->dataHandling == odhTotal) r *= dat->accumulator[bin]->count();
      os << "\t" << r;      

    } else if ( dat->dataType == odtVector3 ) {
      dynamic_cast<VectorAccumulator*>(dat->accumulator[bin])->getStdDev(v);
      if (isinf(v[0]) || isnan(v[0]) || 
          isinf(v[1]) || isnan(v[1]) || 
          isinf(v[2]) || isnan(v[2]) ) {      
        sprintf( painCave.errMsg,
                 "SpatialStatistics detected a numerical error writing:\n"
                 "\tstandard deviation of %s for bin %d",
                 dat->title.c_str(), bin);
        painCave.isFatal = 1;
        simError();
      }
      if (dat->dataHandling == odhTotal) v *= dat->accumulator[bin]->count();
      os << "\t" << v[0] << "\t" << v[1] << "\t" << v[2];
    }
  }
  
  
  OutputData* SpatialStatistics::beginOutputData(vector<OutputData*>::iterator& i) {
    i = data_.begin();
    return i != data_.end()? *i : NULL;
  }

  OutputData* SpatialStatistics::nextOutputData(vector<OutputData*>::iterator& i){
    ++i;
    return i != data_.end()? *i: NULL;
  }


  SlabStatistics::SlabStatistics(SimInfo* info, const string& filename, 
                                 const string& sele, int nbins) : 
    SpatialStatistics(info, filename, sele, nbins) {
    
    z_ = new OutputData;
    z_->units =  "Angstroms";
    z_->title =  "Z";
    z_->dataType = odtReal;
    z_->dataHandling = odhAverage;
    z_->accumulator.reserve(nbins);
    for (int i = 0; i < nbins; i++) 
      z_->accumulator.push_back( new Accumulator() );
    data_.push_back(z_);
  }

  void SlabStatistics::processFrame(int istep) {
    RealType z;
    hmat_ = currentSnapshot_->getHmat();
    for (int i = 0; i < nBins_; i++) {
      z = (((RealType)i + 0.5) / (RealType)nBins_) * hmat_(2,2);
      dynamic_cast<Accumulator*>(z_->accumulator[i])->add(z);
    }
    volume_ = currentSnapshot_->getVolume();

    SpatialStatistics::processFrame(istep);
  }

  int SlabStatistics::getBin(Vector3d pos) {
    currentSnapshot_->wrapVector(pos);
    // which bin is this stuntdouble in?
    // wrapped positions are in the range [-0.5*hmat(2,2), +0.5*hmat(2,2)]
    // Shift molecules by half a box to have bins start at 0
    // The modulo operator is used to wrap the case when we are 
    // beyond the end of the bins back to the beginning.
    return int(nBins_ * (pos.z() / hmat_(2,2) + 0.5)) % nBins_;  
  }


  ShellStatistics::ShellStatistics(SimInfo* info, const string& filename, 
                                   const string& sele, int nbins) : 
    SpatialStatistics(info, filename, sele, nbins){
    
    coordinateOrigin_ = V3Zero;
    binWidth_ = 1.0;

    r_ = new OutputData;
    r_->units =  "Angstroms";
    r_->title =  "R";
    r_->dataType = odtReal;
    r_->dataHandling = odhAverage;
    r_->accumulator.reserve(nbins);
    for (int i = 0; i < nbins; i++)
      r_->accumulator.push_back( new Accumulator() );    
    data_.push_back(r_);

    for (int i = 0; i < nbins; i++) {
      RealType r = (((RealType)i + 0.5) * binWidth_);
      dynamic_cast<Accumulator*>(r_->accumulator[i])->add(r);
    }
  }

  int ShellStatistics::getBin(Vector3d pos) {    
    Vector3d rPos = pos - coordinateOrigin_;
    return int(rPos.length() / binWidth_);
  }
}

