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

#define HONKING_LARGE_VALUE 1.0e10

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

#ifdef IS_MPI
    int nProc, worldRank;

    nProc = MPI::COMM_WORLD.Get_size();
    worldRank = MPI::COMM_WORLD.Get_rank();

    bool my_min_found = min_found;
    bool my_max_found = max_found;

    // Even if we didn't find a minimum, did someone else?
    MPI::COMM_WORLD.Allreduce(&my_min_found, &min_found, 
                              1, MPI::BOOL, MPI::LAND);
    
    // Even if we didn't find a maximum, did someone else?
    MPI::COMM_WORLD.Allreduce(&my_max_found, &max_found, 
                              1, MPI::BOOL, MPI::LAND);
    
    struct {
      RealType val;
      int rank;
    } max_vals, min_vals;
    
    if (min_found) {
      if (my_min_found) 
        min_vals.val = min_val;
      else 
        min_vals.val = HONKING_LARGE_VALUE;
      
      min_vals.rank = worldRank;    
      
      // Who had the minimum?
      MPI::COMM_WORLD.Allreduce(&min_vals, &min_vals, 
                                1, MPI::REALTYPE_INT, MPI::MINLOC);
      min_val = min_vals.val;
    }
      
    if (max_found) {
      if (my_max_found) 
        max_vals.val = max_val;
      else 
        max_vals.val = -HONKING_LARGE_VALUE;
      
      max_vals.rank = worldRank;    
      
      // Who had the maximum?
      MPI::COMM_WORLD.Allreduce(&max_vals, &max_vals, 
                                1, MPI::REALTYPE_INT, MPI::MAXLOC);
      max_val = max_vals.val;
    }
#endif

    if (max_found && min_found) {
      if (min_val< max_val) {

#ifdef IS_MPI       
        if (max_vals.rank == worldRank && min_vals.rank == worldRank) {
          // I have both maximum and minimum, so proceed like a single
          // processor version:
#endif
          // objects to be swapped: velocity & angular velocity
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
#ifdef IS_MPI
          // the rest of the cases only apply in parallel simulations:
        } else if (max_vals.rank == worldRank) {
          // I had the max, but not the minimum
          
          Vector3d min_vel;
          Vector3d max_vel = max_sd->getVel();
          MPI::Status status;

          // point-to-point swap of the velocity vector
          MPI::COMM_WORLD.Sendrecv(max_vel.getArrayPointer(), 3, MPI::REALTYPE,
                                   min_vals.rank, 0, 
                                   min_vel.getArrayPointer(), 3, MPI::REALTYPE,
                                   min_vals.rank, 0, status);
          
          switch(rnemdType_) {
          case rnemdKinetic :
            max_sd->setVel(min_vel);
            
            if (max_sd->isDirectional()) {
              Vector3d min_angMom;
              Vector3d max_angMom = max_sd->getJ();

              // point-to-point swap of the angular momentum vector
              MPI::COMM_WORLD.Sendrecv(max_angMom.getArrayPointer(), 3, 
                                       MPI::REALTYPE, min_vals.rank, 1, 
                                       min_angMom.getArrayPointer(), 3, 
                                       MPI::REALTYPE, min_vals.rank, 1, 
                                       status);

              max_sd->setJ(min_angMom);
            }
            break;
          case rnemdPx :
            max_vel.x() = min_vel.x();
            max_sd->setVel(max_vel);
            break;
          case rnemdPy :
            max_vel.y() = min_vel.y();
            max_sd->setVel(max_vel);
            break;
          case rnemdPz :
            max_vel.z() = min_vel.z();
            max_sd->setVel(max_vel);
            break;
          case rnemdUnknown : 
          default :
            break;
          }
        } else if (min_vals.rank == worldRank) {
          // I had the minimum but not the maximum:
          
          Vector3d max_vel;
          Vector3d min_vel = min_sd->getVel();
          MPI::Status status;
          
          // point-to-point swap of the velocity vector
          MPI::COMM_WORLD.Sendrecv(min_vel.getArrayPointer(), 3, MPI::REALTYPE,
                                   max_vals.rank, 0, 
                                   max_vel.getArrayPointer(), 3, MPI::REALTYPE,
                                   max_vals.rank, 0, status);
          
          switch(rnemdType_) {
          case rnemdKinetic :
            min_sd->setVel(max_vel);
            
            if (min_sd->isDirectional()) {
              Vector3d min_angMom = min_sd->getJ();
              Vector3d max_angMom;

              // point-to-point swap of the angular momentum vector
              MPI::COMM_WORLD.Sendrecv(min_angMom.getArrayPointer(), 3, 
                                       MPI::REALTYPE, max_vals.rank, 1, 
                                       max_angMom.getArrayPointer(), 3, 
                                       MPI::REALTYPE, max_vals.rank, 1, 
                                       status);

              min_sd->setJ(max_angMom);
            }
            break;
          case rnemdPx :
            min_vel.x() = max_vel.x();
            min_sd->setVel(min_vel);
            break;
          case rnemdPy :
            min_vel.y() = max_vel.y();
            min_sd->setVel(min_vel);
            break;
          case rnemdPz :
            min_vel.z() = max_vel.z();
            min_sd->setVel(min_vel);
            break;
          case rnemdUnknown : 
          default :
            break;
          }
        }
#endif
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
    RealType time = currentSnap_->getTime();

    stat[Stats::RNEMD_SWAP_TOTAL] = exchangeSum_;

    seleMan_.setSelectionSet(evaluator_.evaluate());

    int selei;
    StuntDouble* sd;
    int idx;

    std::vector<RealType> valueHist(nBins_, 0.0); // keeps track of what's 
                                                  // being averaged
    std::vector<int> valueCount(nBins_, 0);       // keeps track of the 
                                                  // number of degrees of 
                                                  // freedom being averaged

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
      
      RealType mass = sd->getMass();
      Vector3d vel = sd->getVel();
      RealType value;

      switch(rnemdType_) {
      case rnemdKinetic :
	
	value = mass * (vel[0]*vel[0] + vel[1]*vel[1] + 
			vel[2]*vel[2]);
	
        valueCount[binNo] += 3;
	if (sd->isDirectional()) {
	  Vector3d angMom = sd->getJ();
	  Mat3x3d I = sd->getI();
          
	  if (sd->isLinear()) {
	    int i = sd->linearAxis();
	    int j = (i + 1) % 3;
	    int k = (i + 2) % 3;
	    value += angMom[j] * angMom[j] / I(j, j) + 
	      angMom[k] * angMom[k] / I(k, k);

            valueCount[binNo] +=2;

	  } else {
	    value += angMom[0]*angMom[0]/I(0, 0) 
	      + angMom[1]*angMom[1]/I(1, 1) 
	      + angMom[2]*angMom[2]/I(2, 2);
            valueCount[binNo] +=3;
	  }
	}
	value = value / OOPSEConstant::energyConvert / OOPSEConstant::kb;

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
      valueHist[binNo] += value;
    }

#ifdef IS_MPI

    // all processors have the same number of bins, and STL vectors pack their 
    // arrays, so in theory, this should be safe:

    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &valueHist[0],
                              nBins_, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &valueCount[0],
                              nBins_, MPI::INT, MPI::SUM);

    // If we're the root node, should we print out the results
    int worldRank = MPI::COMM_WORLD.Get_rank();
    if (worldRank == 0) {
#endif
       
      std::cout << time;
      for (int j = 0; j < nBins_; j++)
        std::cout << "\t" << valueHist[j] / (RealType)valueCount[j];
      std::cout << "\n";
      
#ifdef IS_MPI
    }
#endif
  }
}
