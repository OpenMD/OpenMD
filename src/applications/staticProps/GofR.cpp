/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */

#include <algorithm>
#include <fstream>
#include "applications/staticProps/GofR.hpp"
#include "utils/simError.h"

namespace oopse {

GofR::GofR(SimInfo* info, const std::string& filename, const std::string& sele1, const std::string& sele2, double len, int nrbins)
    : RadialDistrFunc(info, filename, sele1, sele2), len_(len), nRBins_(nrbins){

    deltaR_ = len_ /nRBins_;
    
    histogram_.resize(nRBins_);
    avgGofr_.resize(nRBins_);

    setOutputName(getPrefix(filename) + ".gofr");
}


void GofR::preProcess() {
    std::fill(avgGofr_.begin(), avgGofr_.end(), 0.0);    
}

void GofR::initalizeHistogram() {
    std::fill(histogram_.begin(), histogram_.end(), 0);
}


void GofR::processHistogram() {

    int nPairs = getNPairs();
    double volume = info_->getSnapshotManager()->getCurrentSnapshot()->getVolume();
    double pairDensity = nPairs /volume * 2.0;
    double pairConstant = ( 4.0 * NumericConstant::PI * pairDensity ) / 3.0;

    for(int i = 0 ; i < histogram_.size(); ++i){

        double rLower = i * deltaR_;
        double rUpper = rLower + deltaR_;
        double volSlice = ( rUpper * rUpper * rUpper ) - ( rLower * rLower * rLower );
        double nIdeal = volSlice * pairConstant;

        avgGofr_[i] += histogram_[i] / nIdeal;    
    }

}

void GofR::collectHistogram(StuntDouble* sd1, StuntDouble* sd2) {

    if (sd1 == sd2) {
        return;
    }
    
    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12 = pos1 - pos2;
    currentSnapshot_->wrapVector(r12);

    double distance = r12.length();

    if (distance < len_) {
        int whichBin = distance / deltaR_;
        histogram_[whichBin] += 2;
    }
}


void GofR::writeRdf() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
        rdfStream << "#radial distribution function\n";
        rdfStream << "#selection1: (" << selectionScript1_ << ")\t";
        rdfStream << "selection2: (" << selectionScript2_ << ")\n";
        rdfStream << "#r\tcorrValue\n";
        for (int i = 0; i < avgGofr_.size(); ++i) {
            double r = deltaR_ * (i + 0.5);
            rdfStream << r << "\t" << avgGofr_[i]/nProcessed_ << "\n";
        }
        
    } else {


    }

    rdfStream.close();
}

}

