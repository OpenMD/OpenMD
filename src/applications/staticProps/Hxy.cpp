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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [4] , Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). *
 *
 *  Created by Xiuquan Sun on 05/09/06.
 *  @author  Xiuquan Sun 
 *  @version $Id$
 *
 */

/* Calculates the undulation spectrum of the lipid membrance. */

#include <algorithm>
#include <fstream>
#include "applications/staticProps/Hxy.hpp"
#include "utils/simError.h"
#include "utils/PhysicalConstants.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "types/LennardJonesAdapter.hpp"
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

namespace OpenMD {
  
  Hxy::Hxy(SimInfo* info, const std::string& filename,
           const std::string& sele, int nbins_x, int nbins_y, int nbins_z,
           int nrbins)
    : StaticAnalyser(info, filename), selectionScript_(sele),
      evaluator_(info), seleMan_(info), nBinsX_(nbins_x), nBinsY_(nbins_y),
      nBinsZ_(nbins_z), nbins_(nrbins){

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // dens_ stores the local density, rho(x,y,z) on a 3-D grid
    dens_.resize(nBinsX_);
    // bin stores the upper and lower surface cutoff locations (z) for
    // a column through grid location x,y
    minHeight_.resize(nBinsX_);
    maxHeight_.resize(nBinsX_);

    for (unsigned int i = 0; i < nBinsX_; i++) {
      dens_[i].resize(nBinsY_);
      minHeight_[i].resize(nBinsY_);
      maxHeight_[i].resize(nBinsY_);            
      for (unsigned int j = 0; j < nBinsY_; j++) {
        dens_[i][j].resize(nBinsZ_);
      }
    }
    
    mag.resize(nBinsX_*nBinsY_);
    newmag.resize(nBinsX_*nBinsY_);
    sum_bin.resize(nbins_);
    avg_bin.resize(nbins_);
    errbin_sum.resize(nbins_);
    errbin.resize(nbins_);
    sum_bin_sq.resize(nbins_);
    avg_bin_sq.resize(nbins_);
    errbin_sum_sq.resize(nbins_);
    errbin_sq.resize(nbins_);

    bin.resize(nbins_);
    samples.resize(nbins_);

    setOutputName(getPrefix(filename) + ".Hxy");
  }

  Hxy::~Hxy(){
  }

  void Hxy::process() {
#if defined(HAVE_FFTW_H) || defined(HAVE_DFFTW_H) || defined(HAVE_FFTW3_H)
    Molecule* mol;
    RigidBody* rb;
    StuntDouble* sd;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    int ii;

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames/step_;
    for(unsigned int k=0; k < bin.size(); k++)
      bin[k].resize(nFrames);
    for(unsigned int k=0; k < samples.size(); k++)
      samples[k].resize(nFrames);

    for (int istep = 0; istep < nFrames; istep += step_) {

      for (unsigned int i = 0; i < nBinsX_; i++) {
        std::fill(minHeight_[i].begin(), minHeight_[i].end(), 0.0);
        std::fill(maxHeight_[i].begin(), maxHeight_[i].end(), 0.0);
        for (unsigned int j = 0; j < nBinsY_; j++) {
          std::fill(dens_[i][j].begin(), dens_[i][j].end(), 0.0);
        }
      }
                 
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }
      
      // update the positions of atoms which belong to the rigidbodies
      
      for (mol = info_->beginMolecule(mi); mol != NULL; 
           mol = info_->nextMolecule(mi)) {
        for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
             rb = mol->nextRigidBody(rbIter)) {
          rb->updateAtoms();
        }        
      }                       
      
#ifdef HAVE_FFTW3_H
      fftw_plan p;
#else
      fftwnd_plan p;
#endif
      fftw_complex *in, *out;
      
      in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nBinsX_*nBinsY_));
      out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(nBinsX_*nBinsY_));

#ifdef HAVE_FFTW3_H
      p = fftw_plan_dft_2d(nBinsX_, nBinsY_, in, out, 
                           FFTW_FORWARD, FFTW_ESTIMATE); 
#else
      p = fftw2d_create_plan(nBinsX_, nBinsY_, FFTW_FORWARD, FFTW_ESTIMATE);
#endif

      Mat3x3d hmat = currentSnapshot_->getHmat();
      RealType lenX_ = hmat(0,0);
      RealType lenY_ = hmat(1,1);
      RealType lenZ_ = hmat(2,2);
      
      RealType voxelSize = lenX_*lenY_*lenZ_ /
        RealType(nBinsX_*nBinsY_*nBinsZ_);

      RealType x, y, z;
            
      for (sd = seleMan_.beginSelected(ii); sd != NULL;
           sd = seleMan_.nextSelected(ii)) {

        if (sd->isAtom()) {
          Atom* atom = static_cast<Atom*>(sd);
          Vector3d pos = sd->getPos();
          LennardJonesAdapter lja = LennardJonesAdapter(atom->getAtomType());
          RealType mass = sd->getMass();
          RealType sigma = lja.getSigma() * 0.5;
          RealType sigma2 = sigma * sigma;

          for (unsigned int i = 0; i < nBinsX_; i++) {
            x = lenX_ * (RealType(i) / RealType(nBinsX_) );
            
            for (unsigned int j = 0; j < nBinsY_; j++) {
              y = lenY_ * (RealType(j) / RealType(nBinsY_));
              
              for (unsigned int k = 0; k < nBinsZ_; k++) {
                z = lenZ_ * (RealType(k) / RealType(nBinsZ_));
                
                Vector3d r = Vector3d(x, y, z) - pos;
                
                if (usePeriodicBoundaryConditions_)
                  currentSnapshot_->wrapVector(r);
                
                RealType dist = r.length();
                RealType density = mass * exp(-dist*dist/(sigma2*2.0)) /
                  (voxelSize * sqrt(2.0*NumericConstant::PI*sigma2));
                dens_[i][j][k] += PhysicalConstants::densityConvert * density;
              }
            }
          }
        }
      }

      RealType maxDens(0.0);
      for (unsigned int i = 0; i < nBinsX_; i++) {
        for (unsigned int j = 0; j < nBinsY_; j++) {
          for (unsigned int k = 0; k < nBinsZ_; k++) {
            std::cout << dens_[i][j][k] << "\t";
            if (dens_[i][j][k] > maxDens) maxDens = dens_[i][j][k];
          }
          std::cout << "\n";
        }
        std::cout << "\n";
      }
      
      RealType threshold = maxDens / 2.0;

      std::cerr << "threshold = " << threshold << "\n";
      
      for (unsigned int i = 0; i < nBinsX_; i++) {        
        for (unsigned int j = 0; j < nBinsY_; j++) {
          bool lowFound = false;
          for (unsigned int k = 0; k < nBinsZ_; k++) {
            if (!lowFound && dens_[i][j][k] > threshold) {
              minHeight_[i][j] = lenZ_ * (RealType(k) / RealType(nBinsZ_));
              lowFound = true;
            } 
            if (lowFound && dens_[i][j][k] < threshold) {
              maxHeight_[i][j] = lenZ_ * (RealType(k) / RealType(nBinsZ_));
            }                         
          }
        }
      }

     
      int newindex;
      for (unsigned int i=0; i < nBinsX_; i++) {
	for (unsigned int j=0; j < nBinsY_; j++) {
	  newindex = i*nBinsY_ + j;
          //std::cout << maxHeight_[i][j] << "\t";
	  c_re(in[newindex]) = maxHeight_[i][j];
	  c_im(in[newindex]) = 0.0;
	}
        //std::cout << "\n";
      }

#ifdef HAVE_FFTW3_H
      fftw_execute(p);
#else
      fftwnd_one(p, in, out);
#endif
      
      for (unsigned int i=0; i< nBinsX_; i++) {
	for(unsigned int j=0; j< nBinsY_; j++) {
	  newindex = i*nBinsY_ + j;
	  mag[newindex] = pow(c_re(out[newindex]),2) + pow(c_im(out[newindex]),2);
	}
      }

#ifdef HAVE_FFTW3_H
      fftw_destroy_plan(p);
#else
      fftwnd_destroy_plan(p);
#endif      
      fftw_free(out);
      fftw_free(in);

      int index, new_i, new_j, new_index;
      for (unsigned int i=0; i< (nBinsX_/2); i++) {
	for(unsigned int j=0; j< (nBinsY_/2); j++) {
          index = i*nBinsY_ + j;
          new_i = i + (nBinsX_/2);
          new_j = j + (nBinsY_/2);
          new_index = new_i*nBinsY_ + new_j;
	  newmag[new_index] = mag[index];
	}
      }
      
      for (unsigned int i=(nBinsX_/2); i< nBinsX_; i++) {
	for(unsigned int j=0; j< (nBinsY_/2); j++) {
	  index = i*nBinsY_ + j;
	  new_i = i - (nBinsX_/2);
	  new_j = j + (nBinsY_/2);
	  new_index = new_i*nBinsY_ + new_j;
	  newmag[new_index] = mag[index];
	}
      }
      
      for (unsigned int i=0; i< (nBinsX_/2); i++) {
	for(unsigned int j=(nBinsY_/2); j< nBinsY_; j++) {
	  index = i*nBinsY_ + j;
	  new_i = i + (nBinsX_/2);
	  new_j = j - (nBinsY_/2);
	  new_index = new_i*nBinsY_ + new_j;
	  newmag[new_index] = mag[index];
	}
      }
      
      for (unsigned int i=(nBinsX_/2); i< nBinsX_; i++) {
	for(unsigned int j=(nBinsY_/2); j< nBinsY_; j++) {
	  index = i*nBinsY_ + j;
	  new_i = i - (nBinsX_/2);
	  new_j = j - (nBinsY_/2);
	  new_index = new_i*nBinsY_ + new_j;
	  newmag[new_index] = mag[index];
	}
      }
    
      RealType maxfreqx = RealType(nBinsX_) / lenX_;
      RealType maxfreqy = RealType(nBinsY_) / lenY_;
      
      //  printf("%lf\t%lf\t%lf\t%lf\n", dx, dy, maxfreqx, maxfreqy);
      
      RealType maxfreq = sqrt(maxfreqx*maxfreqx + maxfreqy*maxfreqy);
      RealType dfreq = maxfreq/(RealType)(nbins_-1);
    
      //printf("%lf\n", dfreq);
      
      int zero_freq_x = nBinsX_/2; 
      int zero_freq_y = nBinsY_/2; 
      
      for (int i=0; i< nBinsX_; i++) {
	for(int j=0; j< nBinsY_; j++) {
	  
	  RealType freq_x = (RealType)(i - zero_freq_x)*maxfreqx*2 / nBinsX_;
	  RealType freq_y = (RealType)(j - zero_freq_y)*maxfreqy*2 / nBinsY_;
	  
	  RealType freq = sqrt(freq_x*freq_x + freq_y*freq_y);
	  
	  unsigned int whichbin = (unsigned int) (freq / dfreq);
	  newindex = i*nBinsY_ + j;
	  //	printf("%d %d %lf %lf\n", whichbin, newindex, freq, dfreq);
	  bin[whichbin][istep] += newmag[newindex];
	  samples[whichbin][istep]++;
	}
      }
      
      for (unsigned int i = 0; i < nbins_; i++) {
	if ( samples[i][istep] > 0) {
	  bin[i][istep] = 4.0 * sqrt(bin[i][istep] / (RealType)samples[i][istep]) / (RealType)nBinsX_ / (RealType)nBinsY_;
	}
      }    
    }

    for (unsigned int i = 0; i < nbins_; i++) {
      for (unsigned int j = 0; j < nFrames; j++) {
	sum_bin[i] += bin[i][j];
	sum_bin_sq[i] += bin[i][j] * bin[i][j];
      }
      avg_bin[i] = sum_bin[i] / (RealType)nFrames;
      avg_bin_sq[i] = sum_bin_sq[i] / (RealType)nFrames;
      for (int j = 0; j < nFrames; j++) {
	errbin_sum[i] += pow((bin[i][j] - avg_bin[i]), 2);
	errbin_sum_sq[i] += pow((bin[i][j] * bin[i][j] - avg_bin_sq[i]), 2);
      }
      errbin[i] = sqrt( errbin_sum[i] / (RealType)nFrames );
      errbin_sq[i] = sqrt( errbin_sum_sq[i] / (RealType)nFrames );
    }

    printSpectrum();

#else
    sprintf(painCave.errMsg, "Hxy: FFTW support was not compiled in!\n");
    painCave.isFatal = 1;
    simError();  

#endif
  }
    
  void Hxy::printSpectrum() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {

      for (int i = 0; i < nbins_; ++i) {
	if ( avg_bin[i] > 0 ){
	  rdfStream << (RealType)i * dfreq_ << "\t"
                    <<pow(avg_bin[i], 2)<<"\t"
	            <<errbin_sq[i]<<"\t"
	            <<avg_bin[i]<<"\t"
	            <<errbin[i]<<"\n";
	}
      }
    } else {
      
      sprintf(painCave.errMsg, "Hxy: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }

    rdfStream.close();

  }
  
}
