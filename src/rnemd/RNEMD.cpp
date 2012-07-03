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
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */

#include <cmath>
#include "rnemd/RNEMD.hpp"
#include "math/Vector3.hpp"
#include "math/Vector.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Polynomial.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/StuntDouble.hpp"
#include "utils/PhysicalConstants.hpp"
#include "utils/Tuple.hpp"

#ifndef IS_MPI
#include "math/SeqRandNumGen.hpp"
#else
#include "math/ParallelRandNumGen.hpp"
#include <mpi.h>
#endif

#define HONKING_LARGE_VALUE 1.0e10

using namespace std;
namespace OpenMD {
  
  RNEMD::RNEMD(SimInfo* info) : info_(info), evaluator_(info), seleMan_(info), 
                                usePeriodicBoundaryConditions_(info->getSimParams()->getUsePeriodicBoundaryConditions()) {

    failTrialCount_ = 0;
    failRootCount_ = 0;

    int seedValue;
    Globals * simParams = info->getSimParams();
    RNEMDParameters* rnemdParams = simParams->getRNEMDParameters();

    stringToEnumMap_["KineticSwap"] = rnemdKineticSwap;
    stringToEnumMap_["KineticScale"] = rnemdKineticScale;
    stringToEnumMap_["KineticScaleVAM"] = rnemdKineticScaleVAM;
    stringToEnumMap_["KineticScaleAM"] = rnemdKineticScaleAM;
    stringToEnumMap_["PxScale"] = rnemdPxScale;
    stringToEnumMap_["PyScale"] = rnemdPyScale;
    stringToEnumMap_["PzScale"] = rnemdPzScale;
    stringToEnumMap_["Px"] = rnemdPx;
    stringToEnumMap_["Py"] = rnemdPy;
    stringToEnumMap_["Pz"] = rnemdPz;
    stringToEnumMap_["ShiftScaleV"] = rnemdShiftScaleV;
    stringToEnumMap_["ShiftScaleVAM"] = rnemdShiftScaleVAM;
    stringToEnumMap_["Unknown"] = rnemdUnknown;

    runTime_ = simParams->getRunTime();
    statusTime_ = simParams->getStatusTime();

    rnemdObjectSelection_ = rnemdParams->getObjectSelection();
    evaluator_.loadScriptString(rnemdObjectSelection_);
    seleMan_.setSelectionSet(evaluator_.evaluate());

    // do some sanity checking

    int selectionCount = seleMan_.getSelectionCount();
    int nIntegrable = info->getNGlobalIntegrableObjects();

    if (selectionCount > nIntegrable) {
      sprintf(painCave.errMsg, 
              "RNEMD: The current RNEMD_objectSelection,\n"
              "\t\t%s\n"
              "\thas resulted in %d selected objects.  However,\n"
              "\tthe total number of integrable objects in the system\n"
              "\tis only %d.  This is almost certainly not what you want\n"
              "\tto do.  A likely cause of this is forgetting the _RB_0\n"
              "\tselector in the selection script!\n", 
              rnemdObjectSelection_.c_str(), 
              selectionCount, nIntegrable);
      painCave.isFatal = 0;
      painCave.severity = OPENMD_WARNING;
      simError();
    }
    
    const string st = rnemdParams->getExchangeType();

    map<string, RNEMDTypeEnum>::iterator i;
    i = stringToEnumMap_.find(st);
    rnemdType_ = (i == stringToEnumMap_.end()) ? RNEMD::rnemdUnknown : i->second;
    if (rnemdType_ == rnemdUnknown) {
      sprintf(painCave.errMsg, 
              "RNEMD: The current RNEMD_exchangeType,\n"
              "\t\t%s\n"
              "\tis not one of the recognized exchange types.\n",
              st.c_str());
      painCave.isFatal = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }
    
    outputTemp_ = false;
    if (rnemdParams->haveOutputTemperature()) {
      outputTemp_ = rnemdParams->getOutputTemperature();
    } else if ((rnemdType_ == rnemdKineticSwap) ||
	       (rnemdType_ == rnemdKineticScale) ||
	       (rnemdType_ == rnemdKineticScaleVAM) ||
	       (rnemdType_ == rnemdKineticScaleAM)) {
      outputTemp_ = true;
    }
    outputVx_ = false;
    if (rnemdParams->haveOutputVx()) {
      outputVx_ = rnemdParams->getOutputVx();
    } else if ((rnemdType_ == rnemdPx) || (rnemdType_ == rnemdPxScale)) {
      outputVx_ = true;
    }
    outputVy_ = false;
    if (rnemdParams->haveOutputVy()) {
      outputVy_ = rnemdParams->getOutputVy();
    } else if ((rnemdType_ == rnemdPy) || (rnemdType_ == rnemdPyScale)) {
      outputVy_ = true;
    }
    output3DTemp_ = false;
    if (rnemdParams->haveOutputXyzTemperature()) {
      output3DTemp_ = rnemdParams->getOutputXyzTemperature();
    }
    outputRotTemp_ = false;
    if (rnemdParams->haveOutputRotTemperature()) {
      outputRotTemp_ = rnemdParams->getOutputRotTemperature();
    }
    // James put this in.
    outputDen_ = false;
    if (rnemdParams->haveOutputDen()) {
      outputDen_ = rnemdParams->getOutputDen();
    }
    outputAh_ = false;
    if (rnemdParams->haveOutputAh()) {
      outputAh_ = rnemdParams->getOutputAh();
    }    
    outputVz_ = false;
    if (rnemdParams->haveOutputVz()) {
      outputVz_ = rnemdParams->getOutputVz();
    } else if ((rnemdType_ == rnemdPz) || (rnemdType_ == rnemdPzScale)) {
      outputVz_ = true;
    }
    

#ifdef IS_MPI
    if (worldRank == 0) {
#endif

      //may have rnemdWriter separately
      string rnemdFileName;

      if (outputTemp_) {
        rnemdFileName = "temperature.log";
      	tempLog_.open(rnemdFileName.c_str());
      }
      if (outputVx_) {
        rnemdFileName = "velocityX.log";
      	vxzLog_.open(rnemdFileName.c_str());
      }
      if (outputVy_) {
        rnemdFileName = "velocityY.log";
      	vyzLog_.open(rnemdFileName.c_str());
      }

      if (output3DTemp_) {
	rnemdFileName = "temperatureX.log";
	xTempLog_.open(rnemdFileName.c_str());
	rnemdFileName = "temperatureY.log";
	yTempLog_.open(rnemdFileName.c_str());
	rnemdFileName = "temperatureZ.log";
	zTempLog_.open(rnemdFileName.c_str());
      }
      if (outputRotTemp_) {
        rnemdFileName = "temperatureR.log";
      	rotTempLog_.open(rnemdFileName.c_str());
      }
      
      //James put this in
      if (outputDen_) {
	rnemdFileName = "Density.log";
	denLog_.open(rnemdFileName.c_str());
      }
      if (outputAh_) {
	rnemdFileName = "Ah.log";
	AhLog_.open(rnemdFileName.c_str());
      }
      if (outputVz_) {
	rnemdFileName = "velocityZ.log";
	vzzLog_.open(rnemdFileName.c_str());
      }
      logFrameCount_ = 0;
#ifdef IS_MPI
    }
#endif

    set_RNEMD_exchange_time(rnemdParams->getExchangeTime());
    set_RNEMD_nBins(rnemdParams->getNbins());
    midBin_ = nBins_ / 2;
    if (rnemdParams->haveBinShift()) {
      if (rnemdParams->getBinShift()) {
        zShift_ = 0.5 / (RealType)(nBins_);
      } else {
        zShift_ = 0.0;
      }
    } else {
      zShift_ = 0.0;
    }
    //cerr << "I shift slabs by " << zShift_ << " Lz\n";
    //shift slabs by half slab width, maybe useful in heterogeneous systems
    //set to 0.0 if not using it; N/A in status output yet
    if (rnemdParams->haveLogWidth()) {
      set_RNEMD_logWidth(rnemdParams->getLogWidth());
      /*arbitary rnemdLogWidth_, no checking;
      if (rnemdLogWidth_ != nBins_ && rnemdLogWidth_ != midBin_ + 1) {
        cerr << "WARNING! RNEMD_logWidth has abnormal value!\n";
        cerr << "Automaically set back to default.\n";
        rnemdLogWidth_ = nBins_;
      }*/
    } else {
      set_RNEMD_logWidth(nBins_);
    }
    tempHist_.resize(rnemdLogWidth_, 0.0);
    tempCount_.resize(rnemdLogWidth_, 0);
    pxzHist_.resize(rnemdLogWidth_, 0.0);
    //vxzCount_.resize(rnemdLogWidth_, 0);
    pyzHist_.resize(rnemdLogWidth_, 0.0);
    //vyzCount_.resize(rnemdLogWidth_, 0);

    mHist_.resize(rnemdLogWidth_, 0.0);
    xTempHist_.resize(rnemdLogWidth_, 0.0);
    yTempHist_.resize(rnemdLogWidth_, 0.0);
    zTempHist_.resize(rnemdLogWidth_, 0.0);
    xyzTempCount_.resize(rnemdLogWidth_, 0);
    rotTempHist_.resize(rnemdLogWidth_, 0.0);
    rotTempCount_.resize(rnemdLogWidth_, 0);
    // James put this in
    DenHist_.resize(rnemdLogWidth_, 0.0);
    pzzHist_.resize(rnemdLogWidth_, 0.0);

    set_RNEMD_exchange_total(0.0);
    if (rnemdParams->haveTargetFlux()) {
      set_RNEMD_target_flux(rnemdParams->getTargetFlux());
    } else {
      set_RNEMD_target_flux(0.0);
    }
    if (rnemdParams->haveTargetJzKE()) {
      set_RNEMD_target_JzKE(rnemdParams->getTargetJzKE());
    } else {
      set_RNEMD_target_JzKE(0.0);
    }
    if (rnemdParams->haveTargetJzpx()) {
      set_RNEMD_target_jzpx(rnemdParams->getTargetJzpx());
    } else {
      set_RNEMD_target_jzpx(0.0);
    }
    jzp_.x() = targetJzpx_;
    njzp_.x() = -targetJzpx_;
    if (rnemdParams->haveTargetJzpy()) {
      set_RNEMD_target_jzpy(rnemdParams->getTargetJzpy());
    } else {
      set_RNEMD_target_jzpy(0.0);
    }
    jzp_.y() = targetJzpy_;
    njzp_.y() = -targetJzpy_;
    if (rnemdParams->haveTargetJzpz()) {
      set_RNEMD_target_jzpz(rnemdParams->getTargetJzpz());
    } else {
      set_RNEMD_target_jzpz(0.0);
    }
    jzp_.z() = targetJzpz_;
    njzp_.z() = -targetJzpz_;

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
    
#ifdef IS_MPI
    if (worldRank == 0) {
#endif
      
      sprintf(painCave.errMsg, 
              "RNEMD: total failed trials: %d\n",
              failTrialCount_);
      painCave.isFatal = 0;
      painCave.severity = OPENMD_INFO;
      simError();
      
      if (outputTemp_) tempLog_.close();
      if (outputVx_)   vxzLog_.close();
      if (outputVy_)   vyzLog_.close();

      if (rnemdType_ == rnemdKineticScale || rnemdType_ == rnemdPxScale ||
	  rnemdType_ == rnemdPyScale) {
        sprintf(painCave.errMsg, 
                "RNEMD: total root-checking warnings: %d\n",
                failRootCount_);
        painCave.isFatal = 0;
        painCave.severity = OPENMD_INFO;
        simError();
      }
      if (output3DTemp_) {
        xTempLog_.close();
        yTempLog_.close();
        zTempLog_.close();
      }
      if (outputRotTemp_) rotTempLog_.close();
      // James put this in
      if (outputDen_) denLog_.close();
      if (outputAh_)  AhLog_.close();
      if (outputVz_)  vzzLog_.close();
      
#ifdef IS_MPI
    }
#endif
  }

  void RNEMD::doSwap() {

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

      int binNo = int(nBins_ * (pos.z() / hmat(2,2) + zShift_ + 0.5)) % nBins_;


      // if we're in bin 0 or the middleBin
      if (binNo == 0 || binNo == midBin_) {
        
        RealType mass = sd->getMass();
        Vector3d vel = sd->getVel();
        RealType value;

        switch(rnemdType_) {
        case rnemdKineticSwap :
          
          value = mass * vel.lengthSquare();
	  
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
	  } //angular momenta exchange enabled
          //energyConvert temporarily disabled
          //make exchangeSum_ comparable between swap & scale
          //value = value * 0.5 / PhysicalConstants::energyConvert;
          value *= 0.5;
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
	} else { //midBin_
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
    MPI::COMM_WORLD.Allreduce(&my_min_found, &min_found, 1, MPI::BOOL, MPI::LOR);
    // Even if we didn't find a maximum, did someone else?
    MPI::COMM_WORLD.Allreduce(&my_max_found, &max_found, 1, MPI::BOOL, MPI::LOR);
#endif

    if (max_found && min_found) {

#ifdef IS_MPI
      struct {
	RealType val;
	int rank;
      } max_vals, min_vals;
      
      if (my_min_found) {
        min_vals.val = min_val;
      } else {
        min_vals.val = HONKING_LARGE_VALUE;
      }
      min_vals.rank = worldRank;    
      
      // Who had the minimum?
      MPI::COMM_WORLD.Allreduce(&min_vals, &min_vals, 
                                1, MPI::REALTYPE_INT, MPI::MINLOC);
      min_val = min_vals.val;
      
      if (my_max_found) {
        max_vals.val = max_val;
      } else {
        max_vals.val = -HONKING_LARGE_VALUE;
      }
      max_vals.rank = worldRank;    
      
      // Who had the maximum?
      MPI::COMM_WORLD.Allreduce(&max_vals, &max_vals, 
                                1, MPI::REALTYPE_INT, MPI::MAXLOC);
      max_val = max_vals.val;
#endif
      
      if (min_val < max_val) {
	
#ifdef IS_MPI       
        if (max_vals.rank == worldRank && min_vals.rank == worldRank) {
          // I have both maximum and minimum, so proceed like a single
          // processor version:
#endif

          Vector3d min_vel = min_sd->getVel();
          Vector3d max_vel = max_sd->getVel();
          RealType temp_vel;
          
          switch(rnemdType_) {
          case rnemdKineticSwap :
            min_sd->setVel(max_vel);
            max_sd->setVel(min_vel);
	    if (min_sd->isDirectional() && max_sd->isDirectional()) {
              Vector3d min_angMom = min_sd->getJ();
              Vector3d max_angMom = max_sd->getJ();
              min_sd->setJ(max_angMom);
              max_sd->setJ(min_angMom);
	    }//angular momenta exchange enabled
	    //assumes same rigid body identity
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
          case rnemdKineticSwap :
            max_sd->setVel(min_vel);
            //angular momenta exchange enabled
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
          case rnemdKineticSwap :
            min_sd->setVel(max_vel);
            //angular momenta exchange enabled
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
          default :
            break;
          }
        }
#endif
        exchangeSum_ += max_val - min_val;
      } else {        
        sprintf(painCave.errMsg, 
                "RNEMD: exchange NOT performed because min_val > max_val\n");
        painCave.isFatal = 0;
        painCave.severity = OPENMD_INFO;
        simError();        
        failTrialCount_++;
      }
    } else {
      sprintf(painCave.errMsg, 
              "RNEMD: exchange NOT performed because selected object\n"
	      "\tnot present in at least one of the two slabs.\n");
      painCave.isFatal = 0;
      painCave.severity = OPENMD_INFO;
      simError();        
      failTrialCount_++;
    }
    
  }
  
  void RNEMD::doScale() {

    Snapshot* currentSnap_ = info_->getSnapshotManager()->getCurrentSnapshot();
    Mat3x3d hmat = currentSnap_->getHmat();

    seleMan_.setSelectionSet(evaluator_.evaluate());

    int selei;
    StuntDouble* sd;
    int idx;

    vector<StuntDouble*> hotBin, coldBin;

    RealType Phx = 0.0;
    RealType Phy = 0.0;
    RealType Phz = 0.0;
    RealType Khx = 0.0;
    RealType Khy = 0.0;
    RealType Khz = 0.0;
    RealType Khw = 0.0;
    RealType Pcx = 0.0;
    RealType Pcy = 0.0;
    RealType Pcz = 0.0;
    RealType Kcx = 0.0;
    RealType Kcy = 0.0;
    RealType Kcz = 0.0;
    RealType Kcw = 0.0;

    for (sd = seleMan_.beginSelected(selei); sd != NULL; 
         sd = seleMan_.nextSelected(selei)) {

      idx = sd->getLocalIndex();

      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:

      if (usePeriodicBoundaryConditions_)
        currentSnap_->wrapVector(pos);

      // which bin is this stuntdouble in?
      // wrapped positions are in the range [-0.5*hmat(2,2), +0.5*hmat(2,2)]

      int binNo = int(nBins_ * (pos.z() / hmat(2,2) + zShift_ + 0.5)) % nBins_;

      // if we're in bin 0 or the middleBin
      if (binNo == 0 || binNo == midBin_) {
        
        RealType mass = sd->getMass();
        Vector3d vel = sd->getVel();
       
        if (binNo == 0) {
          hotBin.push_back(sd);
          Phx += mass * vel.x();
          Phy += mass * vel.y();
          Phz += mass * vel.z();
          Khx += mass * vel.x() * vel.x();
          Khy += mass * vel.y() * vel.y();
          Khz += mass * vel.z() * vel.z();
	  //if (rnemdType_ == rnemdKineticScaleVAM) {
	  if (sd->isDirectional()) {
	    Vector3d angMom = sd->getJ();
	    Mat3x3d I = sd->getI();
	    if (sd->isLinear()) {
	      int i = sd->linearAxis();
	      int j = (i + 1) % 3;
	      int k = (i + 2) % 3;
	      Khw += angMom[j] * angMom[j] / I(j, j) +
		angMom[k] * angMom[k] / I(k, k);
	    } else {
	      Khw += angMom[0]*angMom[0]/I(0, 0)
		+ angMom[1]*angMom[1]/I(1, 1)
		+ angMom[2]*angMom[2]/I(2, 2);
	    }
	  }
	  //}
        } else { //midBin_
          coldBin.push_back(sd);
          Pcx += mass * vel.x();
          Pcy += mass * vel.y();
          Pcz += mass * vel.z();
          Kcx += mass * vel.x() * vel.x();
          Kcy += mass * vel.y() * vel.y();
          Kcz += mass * vel.z() * vel.z();
	  //if (rnemdType_ == rnemdKineticScaleVAM) {
	  if (sd->isDirectional()) {
	    Vector3d angMom = sd->getJ();
	    Mat3x3d I = sd->getI();
	    if (sd->isLinear()) {
	      int i = sd->linearAxis();
	      int j = (i + 1) % 3;
	      int k = (i + 2) % 3;
	      Kcw += angMom[j] * angMom[j] / I(j, j) +
		angMom[k] * angMom[k] / I(k, k);
	    } else {
	      Kcw += angMom[0]*angMom[0]/I(0, 0)
		+ angMom[1]*angMom[1]/I(1, 1)
		+ angMom[2]*angMom[2]/I(2, 2);
	    }
	  }
	  //}
	}
      }
    }
    
    Khx *= 0.5;
    Khy *= 0.5;
    Khz *= 0.5;
    Khw *= 0.5;
    Kcx *= 0.5;
    Kcy *= 0.5;
    Kcz *= 0.5;
    Kcw *= 0.5;

    // std::cerr << "Khx= " << Khx << "\tKhy= " << Khy << "\tKhz= " << Khz
    // 	      << "\tKhw= " << Khw << "\tKcx= " << Kcx << "\tKcy= " << Kcy
    // 	      << "\tKcz= " << Kcz << "\tKcw= " << Kcw << "\n";
    // std::cerr << "Phx= " << Phx << "\tPhy= " << Phy << "\tPhz= " << Phz
    // 	      << "\tPcx= " << Pcx << "\tPcy= " << Pcy << "\tPcz= " <<Pcz<<"\n";

#ifdef IS_MPI
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Phx, 1, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Phy, 1, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Phz, 1, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Pcx, 1, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Pcy, 1, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Pcz, 1, MPI::REALTYPE, MPI::SUM);

    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Khx, 1, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Khy, 1, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Khz, 1, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Khw, 1, MPI::REALTYPE, MPI::SUM);

    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Kcx, 1, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Kcy, 1, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Kcz, 1, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Kcw, 1, MPI::REALTYPE, MPI::SUM);
#endif

    //solve coldBin coeff's first
    RealType px = Pcx / Phx;
    RealType py = Pcy / Phy;
    RealType pz = Pcz / Phz;
    RealType c, x, y, z;
    bool successfulScale = false;
    if ((rnemdType_ == rnemdKineticScaleVAM) ||
	(rnemdType_ == rnemdKineticScaleAM)) {
      //may need sanity check Khw & Kcw > 0

      if (rnemdType_ == rnemdKineticScaleVAM) {
	c = 1.0 - targetFlux_ / (Kcx + Kcy + Kcz + Kcw);
      } else {
	c = 1.0 - targetFlux_ / Kcw;
      }

      if ((c > 0.81) && (c < 1.21)) {//restrict scaling coefficients
	c = sqrt(c);
	std::cerr << "cold slab scaling coefficient: " << c << endl;
	//now convert to hotBin coefficient
	RealType w = 0.0;
	if (rnemdType_ ==  rnemdKineticScaleVAM) {
	  x = 1.0 + px * (1.0 - c);
	  y = 1.0 + py * (1.0 - c);
	  z = 1.0 + pz * (1.0 - c);
	  /* more complicated way
	     w = 1.0 + (Kcw - Kcw * c * c - (c * c * (Kcx + Kcy + Kcz
	     + Khx * px * px + Khy * py * py + Khz * pz * pz)
	     - 2.0 * c * (Khx * px * (1.0 + px) + Khy * py * (1.0 + py)
	     + Khz * pz * (1.0 + pz)) + Khx * px * (2.0 + px)
	     + Khy * py * (2.0 + py) + Khz * pz * (2.0 + pz)
	     - Kcx - Kcy - Kcz)) / Khw; the following is simpler
	  */
	  if ((fabs(x - 1.0) < 0.1) && (fabs(y - 1.0) < 0.1) &&
	      (fabs(z - 1.0) < 0.1)) {
	    w = 1.0 + (targetFlux_ + Khx * (1.0 - x * x) + Khy * (1.0 - y * y)
		       + Khz * (1.0 - z * z)) / Khw;
	  }//no need to calculate w if x, y or z is out of range
	} else {
	  w = 1.0 + targetFlux_ / Khw;
	}
	if ((w > 0.81) && (w < 1.21)) {//restrict scaling coefficients
	  //if w is in the right range, so should be x, y, z.
	  vector<StuntDouble*>::iterator sdi;
	  Vector3d vel;
	  for (sdi = coldBin.begin(); sdi != coldBin.end(); sdi++) {
	    if (rnemdType_ == rnemdKineticScaleVAM) {
	      vel = (*sdi)->getVel() * c;
	      //vel.x() *= c;
	      //vel.y() *= c;
	      //vel.z() *= c;
	      (*sdi)->setVel(vel);
	    }
	    if ((*sdi)->isDirectional()) {
	      Vector3d angMom = (*sdi)->getJ() * c;
	      //angMom[0] *= c;
	      //angMom[1] *= c;
	      //angMom[2] *= c;
	      (*sdi)->setJ(angMom);
	    }
	  }
	  w = sqrt(w);
	  std::cerr << "xh= " << x << "\tyh= " << y << "\tzh= " << z
		    << "\twh= " << w << endl;
	  for (sdi = hotBin.begin(); sdi != hotBin.end(); sdi++) {
	    if (rnemdType_ == rnemdKineticScaleVAM) {
	      vel = (*sdi)->getVel();
	      vel.x() *= x;
	      vel.y() *= y;
	      vel.z() *= z;
	      (*sdi)->setVel(vel);
	    }
	    if ((*sdi)->isDirectional()) {
	      Vector3d angMom = (*sdi)->getJ() * w;
	      //angMom[0] *= w;
	      //angMom[1] *= w;
	      //angMom[2] *= w;
	      (*sdi)->setJ(angMom);
	    }
	  }
	  successfulScale = true;
	  exchangeSum_ += targetFlux_;
	}
      }
    } else {
      RealType a000, a110, c0, a001, a111, b01, b11, c1;
      switch(rnemdType_) {
      case rnemdKineticScale :
	/* used hotBin coeff's & only scale x & y dimensions
	   RealType px = Phx / Pcx;
	   RealType py = Phy / Pcy;
	   a110 = Khy;
	   c0 = - Khx - Khy - targetFlux_;
	   a000 = Khx;
	   a111 = Kcy * py * py;
	   b11 = -2.0 * Kcy * py * (1.0 + py);
	   c1 = Kcy * py * (2.0 + py) + Kcx * px * ( 2.0 + px) + targetFlux_;
	   b01 = -2.0 * Kcx * px * (1.0 + px);
	   a001 = Kcx * px * px;
	*/
	//scale all three dimensions, let c_x = c_y
	a000 = Kcx + Kcy;
	a110 = Kcz;
	c0 = targetFlux_ - Kcx - Kcy - Kcz;
	a001 = Khx * px * px + Khy * py * py;
	a111 = Khz * pz * pz;
	b01 = -2.0 * (Khx * px * (1.0 + px) + Khy * py * (1.0 + py));
	b11 = -2.0 * Khz * pz * (1.0 + pz);
	c1 = Khx * px * (2.0 + px) + Khy * py * (2.0 + py)
	  + Khz * pz * (2.0 + pz) - targetFlux_;
	break;
      case rnemdPxScale :
	c = 1 - targetFlux_ / Pcx;
	a000 = Kcy;
	a110 = Kcz;
	c0 = Kcx * c * c - Kcx - Kcy - Kcz;
	a001 = py * py * Khy;
	a111 = pz * pz * Khz;
	b01 = -2.0 * Khy * py * (1.0 + py);
	b11 = -2.0 * Khz * pz * (1.0 + pz);
	c1 = Khy * py * (2.0 + py) + Khz * pz * (2.0 + pz)
	  + Khx * (fastpow(c * px - px - 1.0, 2) - 1.0);
	break;
      case rnemdPyScale :
	c = 1 - targetFlux_ / Pcy;
	a000 = Kcx;
	a110 = Kcz;
	c0 = Kcy * c * c - Kcx - Kcy - Kcz;
	a001 = px * px * Khx;
	a111 = pz * pz * Khz;
	b01 = -2.0 * Khx * px * (1.0 + px);
	b11 = -2.0 * Khz * pz * (1.0 + pz);
	c1 = Khx * px * (2.0 + px) + Khz * pz * (2.0 + pz)
	  + Khy * (fastpow(c * py - py - 1.0, 2) - 1.0);
	break;
      case rnemdPzScale ://we don't really do this, do we?
	c = 1 - targetFlux_ / Pcz;
	a000 = Kcx;
	a110 = Kcy;
	c0 = Kcz * c * c - Kcx - Kcy - Kcz;
	a001 = px * px * Khx;
	a111 = py * py * Khy;
	b01 = -2.0 * Khx * px * (1.0 + px);
	b11 = -2.0 * Khy * py * (1.0 + py);
	c1 = Khx * px * (2.0 + px) + Khy * py * (2.0 + py)
	  + Khz * (fastpow(c * pz - pz - 1.0, 2) - 1.0);
	break;
      default :
	break;
      }
      
      RealType v1 = a000 * a111 - a001 * a110;
      RealType v2 = a000 * b01;
      RealType v3 = a000 * b11;
      RealType v4 = a000 * c1 - a001 * c0;
      RealType v8 = a110 * b01;
      RealType v10 = - b01 * c0;
      
      RealType u0 = v2 * v10 - v4 * v4;
      RealType u1 = -2.0 * v3 * v4;
      RealType u2 = -v2 * v8 - v3 * v3 - 2.0 * v1 * v4;
      RealType u3 = -2.0 * v1 * v3;
      RealType u4 = - v1 * v1;
      //rescale coefficients
      RealType maxAbs = fabs(u0);
      if (maxAbs < fabs(u1)) maxAbs = fabs(u1);
      if (maxAbs < fabs(u2)) maxAbs = fabs(u2);
      if (maxAbs < fabs(u3)) maxAbs = fabs(u3);
      if (maxAbs < fabs(u4)) maxAbs = fabs(u4);
      u0 /= maxAbs;
      u1 /= maxAbs;
      u2 /= maxAbs;
      u3 /= maxAbs;
      u4 /= maxAbs;
      //max_element(start, end) is also available.
      Polynomial<RealType> poly; //same as DoublePolynomial poly;
      poly.setCoefficient(4, u4);
      poly.setCoefficient(3, u3);
      poly.setCoefficient(2, u2);
      poly.setCoefficient(1, u1);
      poly.setCoefficient(0, u0);
      vector<RealType> realRoots = poly.FindRealRoots();
      
      vector<RealType>::iterator ri;
      RealType r1, r2, alpha0;
      vector<pair<RealType,RealType> > rps;
      for (ri = realRoots.begin(); ri !=realRoots.end(); ri++) {
	r2 = *ri;
	//check if FindRealRoots() give the right answer
	if ( fabs(u0 + r2 * (u1 + r2 * (u2 + r2 * (u3 + r2 * u4)))) > 1e-6 ) {
	  sprintf(painCave.errMsg, 
		  "RNEMD Warning: polynomial solve seems to have an error!");
	  painCave.isFatal = 0;
	  simError();
	  failRootCount_++;
	}
	//might not be useful w/o rescaling coefficients
	alpha0 = -c0 - a110 * r2 * r2;
	if (alpha0 >= 0.0) {
	  r1 = sqrt(alpha0 / a000);
	  if (fabs(c1 + r1 * (b01 + r1 * a001) + r2 * (b11 + r2 * a111))
	      < 1e-6)
	    { rps.push_back(make_pair(r1, r2)); }
	  if (r1 > 1e-6) { //r1 non-negative
	    r1 = -r1;
	    if (fabs(c1 + r1 * (b01 + r1 * a001) + r2 * (b11 + r2 * a111))
		< 1e-6)
	      { rps.push_back(make_pair(r1, r2)); }
	  }
	}
      }
      // Consider combining together the solving pair part w/ the searching
      // best solution part so that we don't need the pairs vector
      if (!rps.empty()) {
	RealType smallestDiff = HONKING_LARGE_VALUE;
	RealType diff;
	pair<RealType,RealType> bestPair = make_pair(1.0, 1.0);
	vector<pair<RealType,RealType> >::iterator rpi;
	for (rpi = rps.begin(); rpi != rps.end(); rpi++) {
	  r1 = (*rpi).first;
	  r2 = (*rpi).second;
	  switch(rnemdType_) {
	  case rnemdKineticScale :
	    diff = fastpow(1.0 - r1, 2) + fastpow(1.0 - r2, 2)
	      + fastpow(r1 * r1 / r2 / r2 - Kcz/Kcx, 2)
	      + fastpow(r1 * r1 / r2 / r2 - Kcz/Kcy, 2);
	    break;
	  case rnemdPxScale :
	    diff = fastpow(1.0 - r1, 2) + fastpow(1.0 - r2, 2)
	      + fastpow(r1 * r1 / r2 / r2 - Kcz/Kcy, 2);
	    break;
	  case rnemdPyScale :
	    diff = fastpow(1.0 - r1, 2) + fastpow(1.0 - r2, 2)
	      + fastpow(r1 * r1 / r2 / r2 - Kcz/Kcx, 2);
	    break;
	  case rnemdPzScale :
	    diff = fastpow(1.0 - r1, 2) + fastpow(1.0 - r2, 2)
	      + fastpow(r1 * r1 / r2 / r2 - Kcy/Kcx, 2);
	  default :
	    break;
	  }
	  if (diff < smallestDiff) {
	    smallestDiff = diff;
	    bestPair = *rpi;
	  }
	}
#ifdef IS_MPI
	if (worldRank == 0) {
#endif
	  sprintf(painCave.errMsg, 
		  "RNEMD: roots r1= %lf\tr2 = %lf\n",
		  bestPair.first, bestPair.second);
	  painCave.isFatal = 0;
	  painCave.severity = OPENMD_INFO;
	  simError();
#ifdef IS_MPI
	}
#endif
	
	switch(rnemdType_) {
	case rnemdKineticScale :
	  x = bestPair.first;
	  y = bestPair.first;
	  z = bestPair.second;
	  break;
	case rnemdPxScale :
	  x = c;
	  y = bestPair.first;
	  z = bestPair.second;
	  break;
	case rnemdPyScale :
	  x = bestPair.first;
	  y = c;
	  z = bestPair.second;
	  break;
	case rnemdPzScale :
	  x = bestPair.first;
	  y = bestPair.second;
	  z = c;
	  break;          
	default :
	  break;
	}
	vector<StuntDouble*>::iterator sdi;
	Vector3d vel;
	for (sdi = coldBin.begin(); sdi != coldBin.end(); sdi++) {
	  vel = (*sdi)->getVel();
	  vel.x() *= x;
	  vel.y() *= y;
	  vel.z() *= z;
	  (*sdi)->setVel(vel);
	}
	//convert to hotBin coefficient
	x = 1.0 + px * (1.0 - x);
	y = 1.0 + py * (1.0 - y);
	z = 1.0 + pz * (1.0 - z);
	for (sdi = hotBin.begin(); sdi != hotBin.end(); sdi++) {
	  vel = (*sdi)->getVel();
	  vel.x() *= x;
	  vel.y() *= y;
	  vel.z() *= z;
	  (*sdi)->setVel(vel);
	}
	successfulScale = true;
	exchangeSum_ += targetFlux_;
      }
    }
    if (successfulScale != true) {
      sprintf(painCave.errMsg, 
	      "RNEMD: exchange NOT performed!\n");
      painCave.isFatal = 0;
      painCave.severity = OPENMD_INFO;
      simError();        
      failTrialCount_++;
    }
  }

  void RNEMD::doShiftScale() {

    Snapshot* currentSnap_ = info_->getSnapshotManager()->getCurrentSnapshot();
    RealType time = currentSnap_->getTime();	 
    Mat3x3d hmat = currentSnap_->getHmat();

    seleMan_.setSelectionSet(evaluator_.evaluate());

    int selei;
    StuntDouble* sd;
    int idx;

    vector<StuntDouble*> hotBin, coldBin;

    Vector3d Ph(V3Zero);
    RealType Mh = 0.0;
    RealType Kh = 0.0;
    Vector3d Pc(V3Zero);
    RealType Mc = 0.0;
    RealType Kc = 0.0;
    

    for (sd = seleMan_.beginSelected(selei); sd != NULL; 
         sd = seleMan_.nextSelected(selei)) {

      idx = sd->getLocalIndex();

      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:

      if (usePeriodicBoundaryConditions_)
        currentSnap_->wrapVector(pos);

      // which bin is this stuntdouble in?
      // wrapped positions are in the range [-0.5*hmat(2,2), +0.5*hmat(2,2)]

      int binNo = int(nBins_ * (pos.z() / hmat(2,2) + zShift_ + 0.5)) % nBins_;

      // if we're in bin 0 or the middleBin
      if (binNo == 0 || binNo == midBin_) {
        
        RealType mass = sd->getMass();
        Vector3d vel = sd->getVel();
       
        if (binNo == 0) {
          hotBin.push_back(sd);
	  //std::cerr << "before, velocity = " << vel << endl;
	  Ph += mass * vel;
	  //std::cerr << "after, velocity = " << vel << endl;
	  Mh += mass;
	  Kh += mass * vel.lengthSquare();
	  if (rnemdType_ == rnemdShiftScaleVAM) {
	    if (sd->isDirectional()) {
	      Vector3d angMom = sd->getJ();
	      Mat3x3d I = sd->getI();
	      if (sd->isLinear()) {
		int i = sd->linearAxis();
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;
		Kh += angMom[j] * angMom[j] / I(j, j) +
		  angMom[k] * angMom[k] / I(k, k);
	      } else {
		Kh += angMom[0] * angMom[0] / I(0, 0) +
		  angMom[1] * angMom[1] / I(1, 1) +
		  angMom[2] * angMom[2] / I(2, 2);
	      }
	    }
	  }
        } else { //midBin_
          coldBin.push_back(sd);
	  Pc += mass * vel;
	  Mc += mass;
	  Kc += mass * vel.lengthSquare();
	  if (rnemdType_ == rnemdShiftScaleVAM) {
	    if (sd->isDirectional()) {
	      Vector3d angMom = sd->getJ();
	      Mat3x3d I = sd->getI();
	      if (sd->isLinear()) {
		int i = sd->linearAxis();
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;
		Kc += angMom[j] * angMom[j] / I(j, j) +
		  angMom[k] * angMom[k] / I(k, k);
	      } else {
		Kc += angMom[0] * angMom[0] / I(0, 0) +
		  angMom[1] * angMom[1] / I(1, 1) +
		  angMom[2] * angMom[2] / I(2, 2);
	      }
	    }
	  }
	}
      }
    }
    
    Kh *= 0.5;
    Kc *= 0.5;

    // std::cerr << "Mh= " << Mh << "\tKh= " << Kh << "\tMc= " << Mc
    // 	      << "\tKc= " << Kc << endl;
    // std::cerr << "Ph= " << Ph << "\tPc= " << Pc << endl;
    
#ifdef IS_MPI
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Ph[0], 3, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Pc[0], 3, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Mh, 1, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Kh, 1, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Mc, 1, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Kc, 1, MPI::REALTYPE, MPI::SUM);
#endif

    bool successfulExchange = false;
    if ((Mh > 0.0) && (Mc > 0.0)) {//both slabs are not empty
      Vector3d vc = Pc / Mc;
      Vector3d ac = njzp_ / Mc + vc;
      Vector3d acrec = njzp_ / Mc;
      RealType cNumerator = Kc - targetJzKE_ - 0.5 * Mc * ac.lengthSquare();
      if (cNumerator > 0.0) {
	RealType cDenominator = Kc - 0.5 * Mc * vc.lengthSquare();
	if (cDenominator > 0.0) {
	  RealType c = sqrt(cNumerator / cDenominator);
	  if ((c > 0.9) && (c < 1.1)) {//restrict scaling coefficients
	    Vector3d vh = Ph / Mh;
	    Vector3d ah = jzp_ / Mh + vh;
	    Vector3d ahrec = jzp_ / Mh;
	    RealType hNumerator = Kh + targetJzKE_
	      - 0.5 * Mh * ah.lengthSquare();
	    if (hNumerator > 0.0) {
	      RealType hDenominator = Kh - 0.5 * Mh * vh.lengthSquare();
	      if (hDenominator > 0.0) {
		RealType h = sqrt(hNumerator / hDenominator);
		if ((h > 0.9) && (h < 1.1)) {
		  // std::cerr << "cold slab scaling coefficient: " << c << "\n";
		  // std::cerr << "hot slab scaling coefficient: " << h <<  "\n";
		  vector<StuntDouble*>::iterator sdi;
		  Vector3d vel;
		  for (sdi = coldBin.begin(); sdi != coldBin.end(); sdi++) {
		    //vel = (*sdi)->getVel();
		    vel = ((*sdi)->getVel() - vc) * c + ac;
		    (*sdi)->setVel(vel);
		    if (rnemdType_ == rnemdShiftScaleVAM) {
		      if ((*sdi)->isDirectional()) {
			Vector3d angMom = (*sdi)->getJ() * c;
			(*sdi)->setJ(angMom);
		      }
		    }
		  }
		  for (sdi = hotBin.begin(); sdi != hotBin.end(); sdi++) {
		    //vel = (*sdi)->getVel();
		    vel = ((*sdi)->getVel() - vh) * h + ah;
		    (*sdi)->setVel(vel);
		    if (rnemdType_ == rnemdShiftScaleVAM) {
		      if ((*sdi)->isDirectional()) {
			Vector3d angMom = (*sdi)->getJ() * h;
			(*sdi)->setJ(angMom);
		      }
		    }
		  }
		  successfulExchange = true;
		  exchangeSum_ += targetFlux_;
		  // this is a redundant variable for doShiftScale() so that
		  // RNEMD can output one exchange quantity needed in a job.
		  // need a better way to do this.
		  //cerr << "acx =" << ac.x() << "ahx =" << ah.x() << '\n';
		  //cerr << "acy =" << ac.y() << "ahy =" << ah.y() << '\n';
		  //cerr << "acz =" << ac.z() << "ahz =" << ah.z() << '\n';
		  Asum_ += (ahrec.z() - acrec.z());
		  Jsum_ += (jzp_.z()*((1/Mh)+(1/Mc)));
		  AhCount_ = ahrec.z();
		  if (outputAh_) {
		    AhLog_ << time << "   ";
		    AhLog_ << AhCount_;
		    AhLog_ << endl;
		  }		  
		}
	      }
	    }
	  }
	}
      }
    }
    if (successfulExchange != true) {
      //   sprintf(painCave.errMsg, 
      // 	      "RNEMD: exchange NOT performed!\n");
      //   painCave.isFatal = 0;
      //   painCave.severity = OPENMD_INFO;
      //   simError();        
      failTrialCount_++;
    }
  }

  void RNEMD::doRNEMD() {

    switch(rnemdType_) {
    case rnemdKineticScale :
    case rnemdKineticScaleVAM :
    case rnemdKineticScaleAM :
    case rnemdPxScale :
    case rnemdPyScale :
    case rnemdPzScale :
      doScale();
      break;
    case rnemdKineticSwap :
    case rnemdPx :
    case rnemdPy :
    case rnemdPz :
      doSwap();
      break;
    case rnemdShiftScaleV :
    case rnemdShiftScaleVAM :
      doShiftScale();
      break;
    case rnemdUnknown :
    default :
      break;
    }
  }

  void RNEMD::collectData() {

    Snapshot* currentSnap_ = info_->getSnapshotManager()->getCurrentSnapshot();
    Mat3x3d hmat = currentSnap_->getHmat();

    seleMan_.setSelectionSet(evaluator_.evaluate());

    int selei;
    StuntDouble* sd;
    int idx;

    logFrameCount_++;

    // alternative approach, track all molecules instead of only those
    // selected for scaling/swapping:
    /*
    SimInfo::MoleculeIterator miter;
    vector<StuntDouble*>::iterator iiter;
    Molecule* mol;
    StuntDouble* integrableObject;
    for (mol = info_->beginMolecule(miter); mol != NULL;
      mol = info_->nextMolecule(miter))
      integrableObject is essentially sd
        for (integrableObject = mol->beginIntegrableObject(iiter);
             integrableObject != NULL;
             integrableObject = mol->nextIntegrableObject(iiter))
    */
    for (sd = seleMan_.beginSelected(selei); sd != NULL; 
         sd = seleMan_.nextSelected(selei)) {
      
      idx = sd->getLocalIndex();
      
      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:
      
      if (usePeriodicBoundaryConditions_)
        currentSnap_->wrapVector(pos);
      
      // which bin is this stuntdouble in?
      // wrapped positions are in the range [-0.5*hmat(2,2), +0.5*hmat(2,2)]
      
      int binNo = int(rnemdLogWidth_ * (pos.z() / hmat(2,2) + 0.5)) %
	rnemdLogWidth_;
      // no symmetrization allowed due to arbitary rnemdLogWidth_
      /*
      if (rnemdLogWidth_ == midBin_ + 1)
        if (binNo > midBin_)
          binNo = nBins_ - binNo;
      */
      RealType mass = sd->getMass();
      mHist_[binNo] += mass;
      Vector3d vel = sd->getVel();
      RealType value;
      //RealType xVal, yVal, zVal;

      if (outputTemp_) {
	value = mass * vel.lengthSquare();
        tempCount_[binNo] += 3;
	if (sd->isDirectional()) {
	  Vector3d angMom = sd->getJ();
	  Mat3x3d I = sd->getI();
	  if (sd->isLinear()) {
	    int i = sd->linearAxis();
	    int j = (i + 1) % 3;
	    int k = (i + 2) % 3;
	    value += angMom[j] * angMom[j] / I(j, j) + 
	      angMom[k] * angMom[k] / I(k, k);
            tempCount_[binNo] +=2;
	  } else {
	    value += angMom[0] * angMom[0] / I(0, 0) +
	      angMom[1]*angMom[1]/I(1, 1) +
	      angMom[2]*angMom[2]/I(2, 2);
            tempCount_[binNo] +=3;
	  }
	}
	value = value / PhysicalConstants::energyConvert
	  / PhysicalConstants::kb;//may move to getStatus()
	tempHist_[binNo] += value;
      }
      if (outputVx_) {
	value = mass * vel[0];
        //vxzCount_[binNo]++;
	pxzHist_[binNo] += value;
      }
      if (outputVy_) {
	value = mass * vel[1];
        //vyzCount_[binNo]++;
	pyzHist_[binNo] += value;
      }

      if (output3DTemp_) {
	value = mass * vel.x() * vel.x();
	xTempHist_[binNo] += value;
	value = mass * vel.y() * vel.y() / PhysicalConstants::energyConvert
	  / PhysicalConstants::kb;
	yTempHist_[binNo] += value;
	value = mass * vel.z() * vel.z() / PhysicalConstants::energyConvert
	  / PhysicalConstants::kb;
	zTempHist_[binNo] += value;
	xyzTempCount_[binNo]++;
      }
      if (outputRotTemp_) {
	if (sd->isDirectional()) {
	  Vector3d angMom = sd->getJ();
	  Mat3x3d I = sd->getI();
	  if (sd->isLinear()) {
	    int i = sd->linearAxis();
	    int j = (i + 1) % 3;
	    int k = (i + 2) % 3;
	    value = angMom[j] * angMom[j] / I(j, j) + 
	      angMom[k] * angMom[k] / I(k, k);
            rotTempCount_[binNo] +=2;
	  } else {
	    value = angMom[0] * angMom[0] / I(0, 0) +
	      angMom[1] * angMom[1] / I(1, 1) +
	      angMom[2] * angMom[2] / I(2, 2);
            rotTempCount_[binNo] +=3;
	  }
	}
	value = value / PhysicalConstants::energyConvert
	  / PhysicalConstants::kb;//may move to getStatus()
	rotTempHist_[binNo] += value;
      }
      // James put this in.
      if (outputDen_) {
	//value = 1.0;
	DenHist_[binNo] += 1;
      }
      if (outputVz_) {
	value = mass * vel[2];
	//vyzCount_[binNo]++;
	pzzHist_[binNo] += value;
      }     
    }
  }

  void RNEMD::getStarted() {
    collectData();
    /*now can output profile in step 0, but might not be useful;
    Snapshot* currentSnap_ = info_->getSnapshotManager()->getCurrentSnapshot();
    Stats& stat = currentSnap_->statData;
    stat[Stats::RNEMD_EXCHANGE_TOTAL] = exchangeSum_;
    */
    //may output a header for the log file here
    getStatus();
  }

  void RNEMD::getStatus() {

    Snapshot* currentSnap_ = info_->getSnapshotManager()->getCurrentSnapshot();
    RealType time = currentSnap_->getTime();
    //or to be more meaningful, define another item as exchangeSum_ / time
    int j;

#ifdef IS_MPI

    // all processors have the same number of bins, and STL vectors pack their 
    // arrays, so in theory, this should be safe:

    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &mHist_[0],
			      rnemdLogWidth_, MPI::REALTYPE, MPI::SUM);
    if (outputTemp_) {
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &tempHist_[0],
				rnemdLogWidth_, MPI::REALTYPE, MPI::SUM);
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &tempCount_[0],
				rnemdLogWidth_, MPI::INT, MPI::SUM);
    }
    if (outputVx_) {
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &pxzHist_[0],
				rnemdLogWidth_, MPI::REALTYPE, MPI::SUM);
      //MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &vxzCount_[0],
      //			rnemdLogWidth_, MPI::INT, MPI::SUM);
    }
    if (outputVy_) {
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &pyzHist_[0],
				rnemdLogWidth_, MPI::REALTYPE, MPI::SUM);
      //MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &vyzCount_[0],
      //			rnemdLogWidth_, MPI::INT, MPI::SUM);
    }
    if (output3DTemp_) {
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &xTempHist_[0],
                                rnemdLogWidth_, MPI::REALTYPE, MPI::SUM);
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &yTempHist_[0],
                                rnemdLogWidth_, MPI::REALTYPE, MPI::SUM);
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &zTempHist_[0],
                                rnemdLogWidth_, MPI::REALTYPE, MPI::SUM);
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &xyzTempCount_[0],
                                rnemdLogWidth_, MPI::INT, MPI::SUM);
    }
    if (outputRotTemp_) {
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &rotTempHist_[0],
				rnemdLogWidth_, MPI::REALTYPE, MPI::SUM);
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &rotTempCount_[0],
				rnemdLogWidth_, MPI::INT, MPI::SUM);
    }
    // James put this in
    if (outputDen_) {
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &DenHist_[0],
				rnemdLogWidth_, MPI::REALTYPE, MPI::SUM);
    }
    if (outputAh_) {
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &AhCount_,
				1, MPI::REALTYPE, MPI::SUM);
    }
    if (outputVz_) {
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &pzzHist_[0],
				rnemdLogWidth_, MPI::REALTYPE, MPI::SUM);
    }
    
    // If we're the root node, should we print out the results
    int worldRank = MPI::COMM_WORLD.Get_rank();
    if (worldRank == 0) {
#endif

      if (outputTemp_) {
	tempLog_ << time;
	for (j = 0; j < rnemdLogWidth_; j++) {
	  tempLog_ << "\t" << tempHist_[j] / (RealType)tempCount_[j];
	}
	tempLog_ << endl;
      }
      if (outputVx_) {
	vxzLog_ << time;
	for (j = 0; j < rnemdLogWidth_; j++) {
	  vxzLog_ << "\t" << pxzHist_[j] / mHist_[j];
	}
	vxzLog_ << endl;
      }
      if (outputVy_) {
	vyzLog_ << time;
	for (j = 0; j < rnemdLogWidth_; j++) {
	  vyzLog_ << "\t" << pyzHist_[j] / mHist_[j];
	}
	vyzLog_ << endl;
      }

      if (output3DTemp_) {
	RealType temp;
        xTempLog_ << time;
        for (j = 0; j < rnemdLogWidth_; j++) {
	  if (outputVx_)
	    xTempHist_[j] -= pxzHist_[j] * pxzHist_[j] / mHist_[j];
	  temp = xTempHist_[j] / (RealType)xyzTempCount_[j]
	    / PhysicalConstants::energyConvert / PhysicalConstants::kb;
	  xTempLog_ << "\t" << temp;
        }
        xTempLog_ << endl;
        yTempLog_ << time;
        for (j = 0; j < rnemdLogWidth_; j++) {
          yTempLog_ << "\t" << yTempHist_[j] / (RealType)xyzTempCount_[j];
        }
        yTempLog_ << endl;
        zTempLog_ << time;
        for (j = 0; j < rnemdLogWidth_; j++) {
          zTempLog_ << "\t" << zTempHist_[j] / (RealType)xyzTempCount_[j];
        }
        zTempLog_ << endl;
      }
      if (outputRotTemp_) {
	rotTempLog_ << time;
	for (j = 0; j < rnemdLogWidth_; j++) {
	  rotTempLog_ << "\t" << rotTempHist_[j] / (RealType)rotTempCount_[j];
	}
	rotTempLog_ << endl;
      }
      // James put this in.
      Mat3x3d hmat = currentSnap_->getHmat();
      if (outputDen_) {
	denLog_ << time;
	for (j = 0; j < rnemdLogWidth_; j++) {
	  
	  RealType binVol = hmat(0,0) * hmat(1,1) * (hmat(2,2) / float(nBins_));
	  denLog_ << "\t" << DenHist_[j] / (float(logFrameCount_) * binVol);
	}
	denLog_ << endl;
      }
      if (outputVz_) {
	vzzLog_ << time;
	for (j = 0; j < rnemdLogWidth_; j++) {
	  vzzLog_ << "\t" << pzzHist_[j] / mHist_[j];
	}
	vzzLog_ << endl;
      }      
#ifdef IS_MPI
    }
#endif

    for (j = 0; j < rnemdLogWidth_; j++) {
      mHist_[j] = 0.0;
    }
    if (outputTemp_)
      for (j = 0; j < rnemdLogWidth_; j++) {
	tempCount_[j] = 0;
	tempHist_[j] = 0.0;
      }
    if (outputVx_)
      for (j = 0; j < rnemdLogWidth_; j++) {
	//pxzCount_[j] = 0;
	pxzHist_[j] = 0.0;
      }
    if (outputVy_)
      for (j = 0; j < rnemdLogWidth_; j++) {
	//pyzCount_[j] = 0;
	pyzHist_[j] = 0.0;
      }

    if (output3DTemp_)
      for (j = 0; j < rnemdLogWidth_; j++) {
        xTempHist_[j] = 0.0;
        yTempHist_[j] = 0.0;
        zTempHist_[j] = 0.0;
	xyzTempCount_[j] = 0;
      }
    if (outputRotTemp_)
      for (j = 0; j < rnemdLogWidth_; j++) {
	rotTempCount_[j] = 0;
	rotTempHist_[j] = 0.0;
      }
    // James put this in
    if (outputDen_)
      for (j = 0; j < rnemdLogWidth_; j++) {
	//pyzCount_[j] = 0;
	DenHist_[j] = 0.0;
      }
    if (outputVz_)
      for (j = 0; j < rnemdLogWidth_; j++) {
	//pyzCount_[j] = 0;
	pzzHist_[j] = 0.0;
      }    
     // reset the counter
     
    Numcount_++;
    if (Numcount_ > int(runTime_/statusTime_))
      cerr << "time =" << time << "  Asum =" << Asum_ << '\n';
    if (Numcount_ > int(runTime_/statusTime_))
      cerr << "time =" << time << "  Jsum =" << Jsum_ << '\n';
    
    logFrameCount_ = 0;
  }
}

