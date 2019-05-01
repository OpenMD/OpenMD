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
 *  Created by Hemanta Bhattarai on 04/30/19.
 *  @author  Hemanta Bhattarai
 *  @version $Id$
 *
 */

/* Calculates average EAM density profile for selected atom. */

#include <algorithm>
#include <numeric>
#include <fstream>
#include "applications/staticProps/DensityHistogram.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
namespace OpenMD {

  DensityHistogram::DensityHistogram(SimInfo* info, const std::string& filename,
	     const std::string& sele, int nbins)
    : StaticAnalyser(info, filename, nbins), selectionScript_(sele),
      evaluator_(info), seleMan_(info), nBins_(nbins) {

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }



    setOutputName(getPrefix(filename) + ".EAMDensity");
  }

  void DensityHistogram::process() {
    StuntDouble* sd;
    int ii;

    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames/step_;
    vector<vector<RealType>> density(seleMan_.getSelectionCount(),vector<RealType>(nFrames, 0.0 ));


    nProcessed_ = nFrames/step_;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();



      int sele_index = 0;
      for (sd = seleMan_.beginSelected(ii); sd != NULL;
	       sd = seleMan_.nextSelected(ii)) {
           density[sele_index][istep] = sd->getDensity();
           sele_index += 1;
    }



    }
    //finding average density for each selected atom
    vector<RealType> average_density;
    for(int i = 0; i < seleMan_.getSelectionCount(); ++i )
    {
      average_density.push_back(std::accumulate(density[i].begin(),density[i].end(),0.0)/nFrames);

    }
    //sort the average_density
    std::sort (average_density.begin(),average_density.end());


    RealType min = average_density.front();
    RealType max = average_density.back();


    RealType delta_density = (max-min)/(nBins_ - 2 );
    if(delta_density == 0){
      bincenter_.push_back(min);
      histList_.push_back(average_density.size());
    }
    else{
    //fill the center for histogram
    for(int j = 0; j<=nBins_; ++j )
    {
      bincenter_.push_back(min + (j-1) * delta_density);
      histList_.push_back(0);


    }

    //filling up the histogram whith the densities
    int bin_center_pos = 0;
    bool hist_update;
    for(int jj = 0; jj < seleMan_.getSelectionCount(); ++jj){
      hist_update = true;
      while(hist_update){
        if(average_density[jj] >= bincenter_[bin_center_pos] && average_density[jj]< bincenter_[bin_center_pos + 1] ){
          histList_[bin_center_pos] += 1;
          hist_update = false;
        }
        else{
          bin_center_pos++;
          hist_update = true;
        }
      }
    }
  }



/*

    //filling the histogram with densities and count
    int histList_position = 0;
    RealType pre_density = min;
    histList_.push_back(0);
    bincenter_.push_back(min);
    for(int jj = 0; jj < seleMan_.getSelectionCount(); ++jj){
      if(average_density[jj] == pre_density)
        histList_[histList_position] += 1;
      else{
        pre_density = average_density[jj];
        bincenter_.push_back(pre_density);
        histList_.push_back(1);
        histList_position += 1;
      }

    }
*/


    writeDensity();

  }



  void DensityHistogram::writeDensity() {



    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#EAMDensity\n";
      rdfStream << "#nFrames:\t" << nProcessed_ << "\n";
      rdfStream << "#selection: (" << selectionScript_ << ")\n";
      rdfStream << "#" << "Bin_center" << "\tcount\n";
      for (unsigned int i = 0; i < histList_.size(); ++i) {
        rdfStream << bincenter_[i] << "\t"
                  <<  histList_[i]
                  << "\n";
      }

    } else {

      sprintf(painCave.errMsg, "EAMDensity: unable to open %s\n",
	      outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();
  }

}
