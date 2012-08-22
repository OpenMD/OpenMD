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
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

namespace OpenMD {
  
  Hxy::Hxy(SimInfo* info, const std::string& filename, const std::string& sele, int nbins_x, int nbins_y, int nrbins)
    : StaticAnalyser(info, filename), selectionScript_(sele),  evaluator_(info), seleMan_(info), nBinsX_(nbins_x), nBinsY_(nbins_y), nbins_(nrbins){

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    gridsample_.resize(nBinsX_*nBinsY_);
    gridZ_.resize(nBinsX_*nBinsY_);
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
      gridsample_.clear();
      gridZ_.clear();
      sum_bin.clear();
      avg_bin.clear();
      errbin_sum.clear();
      errbin.clear();
      sum_bin_sq.clear();
      avg_bin_sq.clear();
      errbin_sum_sq.clear();
      errbin_sq.clear();
      
      for(unsigned int i=0; i < bin.size(); i++)
	bin[i].clear();
      for(unsigned int i=0; i < samples.size(); i++)
	samples[i].clear();

      mag.clear();
      newmag.clear();
  }

  void Hxy::process() {
#if defined(HAVE_FFTW_H) || defined(HAVE_DFFTW_H) || defined(HAVE_FFTW3_H)
    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames/step_;
    
    for(unsigned int k=0; k < bin.size(); k++)
      bin[k].resize(nFrames);
    for(unsigned int k=0; k < samples.size(); k++)
      samples[k].resize(nFrames);

    RealType lenX_, lenY_;
    RealType gridX_, gridY_;
    RealType halfBoxX_, halfBoxY_;

    int binNoX, binNoY;
    RealType interpsum, value;
    int ninterp, px, py, newp;
    int newx, newy, newindex, index;
    int new_i, new_j, new_index;

    RealType freq_x, freq_y, zero_freq_x, zero_freq_y, freq;
    RealType maxfreqx, maxfreqy, maxfreq;

    int whichbin;
    int nMolecules;

    std::fill(sum_bin.begin(), sum_bin.end(), 0.0);
    std::fill(avg_bin.begin(), avg_bin.end(), 0.0);
    std::fill(errbin_sum.begin(), errbin_sum.end(), 0.0);
    std::fill(errbin.begin(), errbin.end(), 0.0);
    std::fill(sum_bin_sq.begin(), sum_bin_sq.end(), 0.0);
    std::fill(avg_bin_sq.begin(), avg_bin_sq.end(), 0.0);
    std::fill(errbin_sum_sq.begin(), errbin_sum_sq.end(), 0.0);
    std::fill(errbin_sq.begin(), errbin_sq.end(), 0.0);
    
    for(unsigned int i=0; i < bin.size(); i++)
      std::fill(bin[i].begin(), bin[i].end(), 0.0);
    
    for(unsigned int i=0; i < samples.size(); i++)
      std::fill(samples[i].begin(), samples[i].end(), 0);
    
    for (int istep = 0; istep < nFrames; istep += step_) {
      
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      nMolecules = info_->getNGlobalMolecules();
      
      Mat3x3d hmat = currentSnapshot_->getHmat();
      
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

      std::fill(gridsample_.begin(), gridsample_.end(), 0);
      std::fill(gridZ_.begin(), gridZ_.end(), 0.0);
      std::fill(mag.begin(), mag.end(), 0.0);
      std::fill(newmag.begin(), newmag.end(), 0.0);

      int i, j;   
      
      StuntDouble* sd;
      
      lenX_ = hmat(0,0);
      lenY_ = hmat(1,1);
      
      gridX_ = lenX_ /(nBinsX_);
      gridY_ = lenY_ /(nBinsY_);
      
      halfBoxX_ = lenX_ / 2.0;      
      halfBoxY_ = lenY_ / 2.0;      
      
      if (evaluator_.isDynamic()) {
	seleMan_.setSelectionSet(evaluator_.evaluate());
      }
      
      //wrap the stuntdoubles into a cell     
      for (sd = seleMan_.beginSelected(i); sd != NULL; sd = seleMan_.nextSelected(i)) {
	Vector3d pos = sd->getPos();
        if (usePeriodicBoundaryConditions_)
          currentSnapshot_->wrapVector(pos);
	sd->setPos(pos);
      } 
      
      //determine which atom belongs to which grid
      for (sd = seleMan_.beginSelected(i); sd != NULL; sd = seleMan_.nextSelected(i)) {
	Vector3d pos = sd->getPos();
	//int binNo = (pos.z() /deltaR_) - 1;
	int binNoX = (int) ((pos.x() + halfBoxX_) / gridX_);
	int binNoY = (int) ((pos.y() + halfBoxY_) / gridY_);
	//std::cout << "pos.z = " << pos.z() << " halfBoxZ_ = " << halfBoxZ_ << " deltaR_ = "  << deltaR_ << " binNo = " << binNo << "\n";
	gridZ_[binNoX*nBinsY_+binNoY] += pos.z();
	gridsample_[binNoX*nBinsY_+binNoY]++;
      }
      
      // FFT stuff depends on nx and ny, so delay allocation until we have
      // that information
      
      for(i = 0; i < nBinsX_; i++){
	for(j = 0; j < nBinsY_; j++){
	  newindex = i * nBinsY_ + j;
	  if(gridsample_[newindex] > 0){
	    gridZ_[newindex] = gridZ_[newindex] / (RealType)gridsample_[newindex];
	  }
	}
      }
      
      for (i=0; i< nBinsX_; i++) {
	for(j=0; j< nBinsY_; j++) {
	  newindex = i*nBinsY_ + j;
	  if (gridsample_[newindex] == 0) {
	    // interpolate from surrounding points:
	    
	    interpsum = 0.0;
	    ninterp = 0;
	    
	    //point1 = bottom;
	    
	    px = i;
	    py = j - 1;
	    newp = px*nBinsY_ + py;
	    if ((py >= 0) && (gridsample_[newp] > 0)) {
	      interpsum += gridZ_[newp];
	      ninterp++;
	    } 
	    
	    //point2 = top;
	    
	    px = i;
	    py = j + 1;
	    newp = px*nBinsY_ + py;
	    if ((py < nBinsY_) && (gridsample_[newp] > 0)) {
	      interpsum += gridZ_[newp];
	      ninterp++;
	    } 
	    
	    //point3 = left;
	    
	    px = i - 1;
	    py = j;
	    newp = px*nBinsY_ + py;
	    if ((px >= 0) && (gridsample_[newp] > 0)) {
	      interpsum += gridZ_[newp];
	      ninterp++;
	    }
	    
	    //point4 = right;
	    
	    px = i + 1;
	    py = j;
	    newp = px*nBinsY_ + py;
	    if ( (px < nBinsX_ ) && ( gridsample_[newp] > 0 )) { 
              interpsum += gridZ_[newp];
              ninterp++;
            } 
	
	    value = interpsum / (RealType)ninterp;
	    
	    gridZ_[newindex] = value;
	  }
	}
      }
      
      for (i=0; i < nBinsX_; i++) {
	for (j=0; j < nBinsY_; j++) {
	  newindex = i*nBinsY_ + j;
	  
	  c_re(in[newindex]) = gridZ_[newindex];
	  c_im(in[newindex]) = 0.0;
	} 
      }

#ifdef HAVE_FFTW3_H
      fftw_execute(p);
#else
      fftwnd_one(p, in, out);
#endif
      
      for (i=0; i< nBinsX_; i++) {
	for(j=0; j< nBinsY_; j++) {
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

      for (i=0; i< (nBinsX_/2); i++) {
	for(j=0; j< (nBinsY_/2); j++) {
	  index = i*nBinsY_ + j;
	  new_i = i + (nBinsX_/2);
	  new_j = j + (nBinsY_/2);
	  new_index = new_i*nBinsY_ + new_j;
	  newmag[new_index] = mag[index];
	}
      }
      
      for (i=(nBinsX_/2); i< nBinsX_; i++) {
	for(j=0; j< (nBinsY_/2); j++) {
	  index = i*nBinsY_ + j;
	  new_i = i - (nBinsX_/2);
	  new_j = j + (nBinsY_/2);
	  new_index = new_i*nBinsY_ + new_j;
	  newmag[new_index] = mag[index];
	}
      }
      
      for (i=0; i< (nBinsX_/2); i++) {
	for(j=(nBinsY_/2); j< nBinsY_; j++) {
	  index = i*nBinsY_ + j;
	  new_i = i + (nBinsX_/2);
	  new_j = j - (nBinsY_/2);
	  new_index = new_i*nBinsY_ + new_j;
	  newmag[new_index] = mag[index];
	}
      }
      
      for (i=(nBinsX_/2); i< nBinsX_; i++) {
	for(j=(nBinsY_/2); j< nBinsY_; j++) {
	  index = i*nBinsY_ + j;
	  new_i = i - (nBinsX_/2);
	  new_j = j - (nBinsY_/2);
	  new_index = new_i*nBinsY_ + new_j;
	  newmag[new_index] = mag[index];
	}
      }
    
      maxfreqx = 1.0 / gridX_;
      maxfreqy = 1.0 / gridY_;
      
      //  printf("%lf\t%lf\t%lf\t%lf\n", dx, dy, maxfreqx, maxfreqy);
      
      maxfreq = sqrt(maxfreqx*maxfreqx + maxfreqy*maxfreqy);
      dfreq = maxfreq/(RealType)(nbins_-1);
    
      //printf("%lf\n", dfreq);
      
      zero_freq_x = nBinsX_/2; 
      zero_freq_y = nBinsY_/2; 
      
      for (i=0; i< nBinsX_; i++) {
	for(j=0; j< nBinsY_; j++) {
	  
	  freq_x = (RealType)(i - zero_freq_x)*maxfreqx*2 / nBinsX_;
	  freq_y = (RealType)(j - zero_freq_y)*maxfreqy*2 / nBinsY_;
	  
	  freq = sqrt(freq_x*freq_x + freq_y*freq_y);
	  
	  whichbin = (int) (freq / dfreq);
	  newindex = i*nBinsY_ + j;
	  //	printf("%d %d %lf %lf\n", whichbin, newindex, freq, dfreq);
	  bin[whichbin][istep] += newmag[newindex];
	  samples[whichbin][istep]++;
	}
      }
      
      for ( i = 0; i < nbins_; i++) {
	if ( samples[i][istep] > 0) {
	  bin[i][istep] = 4.0 * sqrt(bin[i][istep] / (RealType)samples[i][istep]) / (RealType)nBinsX_ / (RealType)nBinsY_;
	}
      }    
    }

    for (int i = 0; i < nbins_; i++) {
      for (int j = 0; j < nFrames; j++) {
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
	  rdfStream << (RealType)i * dfreq << "\t"
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
