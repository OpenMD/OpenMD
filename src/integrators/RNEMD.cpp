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
#include "math/Vector3.hpp"
#include "math/SquareMatrix3.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/StuntDouble.hpp"
#include "utils/OOPSEConstant.hpp"
#include "utils/Tuple.hpp"

#ifndef IS_MPI
#include "math/SeqRandNumGen.hpp"
#else
#include "math/ParallelRandNumGen.hpp"
#endif


namespace oopse {
  
  RNEMD::RNEMD(SimInfo* info) : info_(info), evaluator_(info), seleMan_(info), usePeriodicBoundaryConditions_(info->getSimParams()->getUsePeriodicBoundaryConditions()) {
    
    int seedValue;
    Globals * simParams = info->getSimParams();

    stringToEnumMap_["Kinetic"] = rnemdKinetic;
    stringToEnumMap_["Px"] = rnemdPx;
    stringToEnumMap_["Py"] = rnemdPy;
    stringToEnumMap_["Pz"] = rnemdPz;
    stringToEnumMap_["Unknown"] = rnemdUnknown;

    rnemdObjectSelection_ = simParams->getRNEMD_objectSelection();
    evaluator_.loadScriptString(rnemdObjectSelection_);
    seleMan_.setSelectionSet(evaluator_.evaluate());


    // do some sanity checking

    int selectionCount = seleMan_.getSelectionCount();
    int nIntegrable = info->getNGlobalIntegrableObjects();

    if (selectionCount > nIntegrable) {
      sprintf(painCave.errMsg, 
              "RNEMD warning: The current RNEMD_objectSelection,\n"
              "\t\t%s\n"
              "\thas resulted in %d selected objects.  However,\n"
              "\tthe total number of integrable objects in the system\n"
              "\tis only %d.  This is almost certainly not what you want\n"
              "\tto do.  A likely cause of this is forgetting the _RB_0\n"
              "\tselector in the selection script!\n", 
              rnemdObjectSelection_.c_str(), 
              selectionCount, nIntegrable);
      painCave.isFatal = 0;
      simError();

    }
    
    const std::string st = simParams->getRNEMD_swapType();

    std::map<std::string, RNEMDTypeEnum>::iterator i;
    i = stringToEnumMap_.find(st);
    rnemdType_  = (i == stringToEnumMap_.end()) ? RNEMD::rnemdUnknown : i->second;

    set_RNEMD_swapTime(simParams->getRNEMD_swapTime());
    set_RNEMD_nBins(simParams->getRNEMD_nBins());
    exchangeSum_ = 0.0;
    counter_ = 0; //added by shenyu

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
    int midBin = nBins_ / 2;

    Snapshot* currentSnap_ = info_->getSnapshotManager()->getCurrentSnapshot();
    Mat3x3d hmat = currentSnap_->getHmat();

    seleMan_.setSelectionSet(evaluator_.evaluate());

    int selei;
    StuntDouble* sd;
    int idx;

    RealType min_val;
    bool min_found = false;   
    StuntDouble* min_sd;

    RealType max_val;
    bool max_found = false;
    StuntDouble* max_sd;

    for (sd = seleMan_.beginSelected(selei); sd != NULL; 
         sd = seleMan_.nextSelected(selei)) {

      idx = sd->getLocalIndex();

      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:

      if (usePeriodicBoundaryConditions_)
        currentSnap_->wrapVector(pos);

      // which bin is this stuntdouble in?
      // wrapped positions are in the range [-0.5*hmat(2,2), +0.5*hmat(2,2)]

      int binNo = int(nBins_ * (pos.z() / hmat(2,2) + 0.5)) % nBins_;


      // if we're in bin 0 or the middleBin
      if (binNo == 0 || binNo == midBin) {
        
        RealType mass = sd->getMass();
        Vector3d vel = sd->getVel();
        RealType value;

        switch(rnemdType_) {
        case rnemdKinetic :
          
          value = mass * (vel[0]*vel[0] + vel[1]*vel[1] + 
                          vel[2]*vel[2]);
          if (sd->isDirectional()) {
            Vector3d angMom = sd->getJ();
            Mat3x3d I = sd->getI();
            
            if (sd->isLinear()) {
              int i = sd->linearAxis();
              int j = (i + 1) % 3;
              int k = (i + 2) % 3;
              value += angMom[j] * angMom[j] / I(j, j) + 
                angMom[k] * angMom[k] / I(k, k);
            } else {                        
              value += angMom[0]*angMom[0]/I(0, 0) 
                + angMom[1]*angMom[1]/I(1, 1) 
                + angMom[2]*angMom[2]/I(2, 2);
            }
          }
          value = value * 0.5 / OOPSEConstant::energyConvert;
          break;
        case rnemdPx :
          value = mass * vel[0];
          break;
        case rnemdPy :
          value = mass * vel[1];
          break;
        case rnemdPz :
          value = mass * vel[2];
          break;
        case rnemdUnknown : 
        default :
          break;
        }
        
        if (binNo == 0) {
	  if (!min_found) {
	    min_val = value;
	    min_sd = sd;
	    min_found = true;
	  } else {
	    if (min_val > value) {
	      min_val = value;
	      min_sd = sd;
	    }
	  }
	} else {
	  if (!max_found) {
	    max_val = value;
	    max_sd = sd;
	    max_found = true;
	  } else {
	    if (max_val < value) {
	      max_val = value;
	      max_sd = sd;
	    }
	  }	  
	}
      }
    }

    // missing:  swap information in parallel

    if (max_found && min_found) {
      if (min_val< max_val) {

	Vector3d min_vel = min_sd->getVel();
	Vector3d max_vel = max_sd->getVel();
	RealType temp_vel;

	switch(rnemdType_) {
	case rnemdKinetic :
	  min_sd->setVel(max_vel);
	  max_sd->setVel(min_vel);

	  if (min_sd->isDirectional() && max_sd->isDirectional()) {
	    Vector3d min_angMom = min_sd->getJ();
	    Vector3d max_angMom = max_sd->getJ();
	    min_sd->setJ(max_angMom);
	    max_sd->setJ(min_angMom);
          }
	  break;
	case rnemdPx :
	  temp_vel = min_vel.x();
	  min_vel.x() = max_vel.x();
	  max_vel.x() = temp_vel;
	  min_sd->setVel(min_vel);
	  max_sd->setVel(max_vel);
	  break;
	case rnemdPy :
	  temp_vel = min_vel.y();
	  min_vel.y() = max_vel.y();
	  max_vel.y() = temp_vel;
	  min_sd->setVel(min_vel);
	  max_sd->setVel(max_vel);
	  break;
	case rnemdPz :
	  temp_vel = min_vel.z();
	  min_vel.z() = max_vel.z();
	  max_vel.z() = temp_vel;
	  min_sd->setVel(min_vel);
	  max_sd->setVel(max_vel);
	  break;
	case rnemdUnknown : 
	default :
	  break;
	}
      exchangeSum_ += max_val - min_val;
      } else {
	std::cerr << "exchange NOT performed.\nmin_val > max_val.\n";
      }
    } else {
      std::cerr << "exchange NOT performed.\none of the two slabs empty.\n";
    }
  }

  void RNEMD::getStatus() {

    Snapshot* currentSnap_ = info_->getSnapshotManager()->getCurrentSnapshot();
    Mat3x3d hmat = currentSnap_->getHmat();
    Stats& stat = currentSnap_->statData;

    stat[Stats::RNEMD_SWAP_TOTAL] = exchangeSum_;

    seleMan_.setSelectionSet(evaluator_.evaluate());

    int selei;
    StuntDouble* sd;
    int idx;

    std::vector<RealType> valueHist;  // keeps track of what's being averaged
    std::vector<int> valueCount; // keeps track of the number of degrees of 
                                 // freedom being averaged
    valueHist.resize(nBins_);
    valueCount.resize(nBins_);
    //do they initialize themselves to zero automatically?
    for (sd = seleMan_.beginSelected(selei); sd != NULL; 
         sd = seleMan_.nextSelected(selei)) {
      
      idx = sd->getLocalIndex();
      
      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:
      
      if (usePeriodicBoundaryConditions_)
        currentSnap_->wrapVector(pos);
      
      // which bin is this stuntdouble in?
      // wrapped positions are in the range [-0.5*hmat(2,2), +0.5*hmat(2,2)]
      
      int binNo = int(nBins_ * (pos.z() / hmat(2,2) + 0.5)) % nBins_;
      
      //std::cerr << "pos.z() = " << pos.z() << " bin = " << binNo << "\n";
      
      RealType mass = sd->getMass();
      Vector3d vel = sd->getVel();
      //std::cerr << "mass = " << mass << " vel = " << vel << "\n";
      RealType value;

      switch(rnemdType_) {
      case rnemdKinetic :
	
	value = mass * (vel[0]*vel[0] + vel[1]*vel[1] + 
			vel[2]*vel[2]);
	
        valueCount[binNo] += 3;
	//std::cerr <<"starting value = " << value << "\n";
	if (sd->isDirectional()) {
          //std::cerr << "angMom calculated.\n";
	  Vector3d angMom = sd->getJ();
          //std::cerr << "current angMom: " << angMom << "\n";
	  Mat3x3d I = sd->getI();
          
	  if (sd->isLinear()) {
	    int i = sd->linearAxis();
	    int j = (i + 1) % 3;
	    int k = (i + 2) % 3;
	    value += angMom[j] * angMom[j] / I(j, j) + 
	      angMom[k] * angMom[k] / I(k, k);

            valueCount[binNo] +=2;

	  } else {
            //std::cerr << "non-linear molecule.\n";
	    value += angMom[0]*angMom[0]/I(0, 0) 
	      + angMom[1]*angMom[1]/I(1, 1) 
	      + angMom[2]*angMom[2]/I(2, 2);
            valueCount[binNo] +=3;

	  }
	}
	//std::cerr <<"total value = " << value << "\n";
	//value *= 0.5 / OOPSEConstant::energyConvert;  // get it in kcal / mol
        //value *= 2.0 / OOPSEConstant::kb;            // convert to temperature
	value = value / OOPSEConstant::energyConvert / OOPSEConstant::kb;
        //std::cerr <<"value = " << value << "\n";
	break;
      case rnemdPx :
	value = mass * vel[0];
        valueCount[binNo]++;
	break;
      case rnemdPy :
	value = mass * vel[1];
        valueCount[binNo]++;
	break;
      case rnemdPz :
	value = mass * vel[2];
        valueCount[binNo]++;
	break;
      case rnemdUnknown : 
      default :
	break;
      }
      //std::cerr << "bin = " << binNo << " value = " << value ;
      valueHist[binNo] += value;
      //std::cerr << " hist = " << valueHist[binNo] << " count = " << valueCount[binNo] << "\n";
    }
    
    std::cout << counter_++;
    for (int j = 0; j < nBins_; j++)
      std::cout << "\t" << valueHist[j] / (RealType)valueCount[j];
    std::cout << "\n";
  }
}
