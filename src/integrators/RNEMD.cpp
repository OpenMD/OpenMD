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

#include "integrators/RNEMD.hpp"
#include "math/SquareMatrix3.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/StuntDouble.hpp"

#ifndef IS_MPI
#include "math/SeqRandNumGen.hpp"
#else
#include "math/ParallelRandNumGen.hpp"
#endif

/* Remove me after testing*/
/*
#include <cstdio>
#include <iostream>
*/
/*End remove me*/

namespace oopse {
  
  RNEMD::RNEMD(SimInfo* info) : info_(info), evaluator_(info), seleMan_(info) {
    
    int seedValue;
    Globals * simParams = info->getSimParams();

    stringToEnumMap_["Kinetic"] = rnemdKinetic;
    stringToEnumMap_["Px"] = rnemdPx;
    stringToEnumMap_["Py"] = rnemdPy;
    stringToEnumMap_["Pz"] = rnemdPz;
    stringToEnumMap_["Unknown"] = rnemdUnknown;

    rnemdObjectSelection_ = simParams->getRNEMD_objectSelection();

    std::cerr << "calling  evaluator with " << rnemdObjectSelection_ << "\n";
    evaluator_.loadScriptString(rnemdObjectSelection_);
    std::cerr << "done";
    
    const std::string st = simParams->getRNEMD_swapType();

    std::map<std::string, RNEMDTypeEnum>::iterator i;
    i = stringToEnumMap_.find(st);
    rnemdType_  = (i == stringToEnumMap_.end()) ? RNEMD::rnemdUnknown : i->second;

    set_RNEMD_swapTime(simParams->getRNEMD_swapTime());
    set_RNEMD_nBins(simParams->getRNEMD_nBins());
    exchangeSum_ = 0.0;
    
#ifndef IS_MPI
    if (simParams->haveSeed()) {
      seedValue = simParams->getSeed();
      randNumGen_ = new SeqRandNumGen(seedValue);
    }else {
      randNumGen_ = new SeqRandNumGen();
    }    
#else
    if (simParams->haveSeed()) {
      seedValue = simParams->getSeed();
      randNumGen_ = new ParallelRandNumGen(seedValue);
    }else {
      randNumGen_ = new ParallelRandNumGen();
    }    
#endif 
  }
  
  RNEMD::~RNEMD() {
    delete randNumGen_;
  }

  void RNEMD::doSwap() {
    std::cerr << "in RNEMD!\n";   
    std::cerr << "nBins = " << nBins_ << "\n";
    std::cerr << "swapTime = " << swapTime_ << "\n";
    std::cerr << "exchangeSum = " << exchangeSum_ << "\n";
    std::cerr << "swapType = " << rnemdType_ << "\n";
    std::cerr << "selection = " << rnemdObjectSelection_ << "\n";

    seleMan_.setSelectionSet(evaluator_.evaluate());

    std::cerr << "selectionCount = " << seleMan_.getSelectionCount() << "\n\n";

    int i;
    StuntDouble* sd;

    for (sd = seleMan_.beginSelected(i); sd != NULL; sd = seleMan_.nextSelected(i)) {
      Vector3d pos = sd->getPos();
      //wrap the stuntdoubles into a cell      
      if (usePeriodicBoundaryConditions_)
        info_->getSnapshotManager()->getCurrentSnapshot()->wrapVector(pos);
      int binNo = int(nBins_ * (pos.z()) / hmat(2,2));
      sliceSDLists_[binNo].push_back(sd);
    }

  }   
}
