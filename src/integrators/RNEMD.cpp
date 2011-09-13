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
#include "integrators/RNEMD.hpp"
#include "math/Vector3.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Polynomial.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/StuntDouble.hpp"
#include "utils/PhysicalConstants.hpp"
#include "utils/Tuple.hpp"

#ifndef IS_MPI
#include "math/SeqRandNumGen.hpp"
#else
#include <mpi.h>
#include "math/ParallelRandNumGen.hpp"
#endif

#define HONKING_LARGE_VALUE 1.0e10

namespace OpenMD {
  
  RNEMD::RNEMD(SimInfo* info) : info_(info), evaluator_(info), seleMan_(info), usePeriodicBoundaryConditions_(info->getSimParams()->getUsePeriodicBoundaryConditions()) {

    failTrialCount_ = 0;
    failRootCount_ = 0;

    int seedValue;
    Globals * simParams = info->getSimParams();

    stringToEnumMap_["KineticSwap"] = rnemdKineticSwap;
    stringToEnumMap_["KineticScale"] = rnemdKineticScale;
    stringToEnumMap_["PxScale"] = rnemdPxScale;
    stringToEnumMap_["PyScale"] = rnemdPyScale;
    stringToEnumMap_["PzScale"] = rnemdPzScale;
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
    
    const std::string st = simParams->getRNEMD_exchangeType();

    std::map<std::string, RNEMDTypeEnum>::iterator i;
    i = stringToEnumMap_.find(st);
    rnemdType_ = (i == stringToEnumMap_.end()) ? RNEMD::rnemdUnknown : i->second;
    if (rnemdType_ == rnemdUnknown) {
      std::cerr << "WARNING! RNEMD Type Unknown!\n";
    }

#ifdef IS_MPI
    if (worldRank == 0) {
#endif

      std::string rnemdFileName;
      std::string xTempFileName;
      std::string yTempFileName;
      std::string zTempFileName;
      switch(rnemdType_) {
      case rnemdKineticSwap :
      case rnemdKineticScale :
        rnemdFileName = "temperature.log";
        break;
      case rnemdPx :
      case rnemdPxScale :
      case rnemdPy :
      case rnemdPyScale :
        rnemdFileName = "momemtum.log";
        xTempFileName = "temperatureX.log";
        yTempFileName = "temperatureY.log";
        zTempFileName = "temperatureZ.log";
        xTempLog_.open(xTempFileName.c_str());
        yTempLog_.open(yTempFileName.c_str());
        zTempLog_.open(zTempFileName.c_str());
        break;
      case rnemdPz :
      case rnemdPzScale :
      case rnemdUnknown :
      default :
        rnemdFileName = "rnemd.log";
        break;
      }
      rnemdLog_.open(rnemdFileName.c_str());

#ifdef IS_MPI
    }
#endif

    set_RNEMD_exchange_time(simParams->getRNEMD_exchangeTime());
    set_RNEMD_nBins(simParams->getRNEMD_nBins());
    midBin_ = nBins_ / 2;
    if (simParams->haveRNEMD_logWidth()) {
      rnemdLogWidth_ = simParams->getRNEMD_logWidth();
      if (rnemdLogWidth_ != nBins_ && rnemdLogWidth_ != midBin_ + 1) {
        std::cerr << "WARNING! RNEMD_logWidth has abnormal value!\n";
        std::cerr << "Automaically set back to default.\n";
        rnemdLogWidth_ = nBins_;
      }
    } else {
      rnemdLogWidth_ = nBins_;
    }
    valueHist_.resize(rnemdLogWidth_, 0.0);
    valueCount_.resize(rnemdLogWidth_, 0);
    xTempHist_.resize(rnemdLogWidth_, 0.0);
    yTempHist_.resize(rnemdLogWidth_, 0.0);
    zTempHist_.resize(rnemdLogWidth_, 0.0);

    set_RNEMD_exchange_total(0.0);
    if (simParams->haveRNEMD_targetFlux()) {
      set_RNEMD_target_flux(simParams->getRNEMD_targetFlux());
    } else {
      set_RNEMD_target_flux(0.0);
    }

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
      std::cerr << "total fail trials: " << failTrialCount_ << "\n";
      rnemdLog_.close();
      if (rnemdType_ == rnemdKineticScale || rnemdType_ == rnemdPxScale || rnemdType_ == rnemdPyScale)
        std::cerr<< "total root-checking warnings: " << failRootCount_ << "\n";
      if (rnemdType_ == rnemdPx || rnemdType_ == rnemdPxScale || rnemdType_ == rnemdPy || rnemdType_ == rnemdPyScale) {
        xTempLog_.close();
        yTempLog_.close();
        zTempLog_.close();
      }
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

      int binNo = int(nBins_ * (pos.z() / hmat(2,2) + 0.5)) % nBins_;


      // if we're in bin 0 or the middleBin
      if (binNo == 0 || binNo == midBin_) {
        
        RealType mass = sd->getMass();
        Vector3d vel = sd->getVel();
        RealType value;

        switch(rnemdType_) {
        case rnemdKineticSwap :
          
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
          //make exchangeSum_ comparable between swap & scale
          //temporarily without using energyConvert
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
          case rnemdKineticSwap :
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
        std::cerr << "exchange NOT performed!\nmin_val > max_val.\n";
        failTrialCount_++;
      }
    } else {
      std::cerr << "exchange NOT performed!\n";
      std::cerr << "at least one of the two slabs empty.\n";
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

    std::vector<StuntDouble*> hotBin, coldBin;

    RealType Phx = 0.0;
    RealType Phy = 0.0;
    RealType Phz = 0.0;
    RealType Khx = 0.0;
    RealType Khy = 0.0;
    RealType Khz = 0.0;
    RealType Pcx = 0.0;
    RealType Pcy = 0.0;
    RealType Pcz = 0.0;
    RealType Kcx = 0.0;
    RealType Kcy = 0.0;
    RealType Kcz = 0.0;

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
        } else { //midBin_
          coldBin.push_back(sd);
          Pcx += mass * vel.x();
          Pcy += mass * vel.y();
          Pcz += mass * vel.z();
          Kcx += mass * vel.x() * vel.x();
          Kcy += mass * vel.y() * vel.y();
          Kcz += mass * vel.z() * vel.z();
	}
      }
    }

    Khx *= 0.5;
    Khy *= 0.5;
    Khz *= 0.5;
    Kcx *= 0.5;
    Kcy *= 0.5;
    Kcz *= 0.5;

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
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Kcx, 1, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Kcy, 1, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &Kcz, 1, MPI::REALTYPE, MPI::SUM);
#endif

    //use coldBin coeff's
    RealType px = Pcx / Phx;
    RealType py = Pcy / Phy;
    RealType pz = Pcz / Phz;

    RealType a000, a110, c0, a001, a111, b01, b11, c1, c;
    switch(rnemdType_) {
    case rnemdKineticScale :
    /*used hotBin coeff's & only scale x & y dimensions
      RealType px = Phx / Pcx;
      RealType py = Phy / Pcy;
      a110 = Khy;
      c0 = - Khx - Khy - targetFlux_;
      a000 = Khx;
      a111 = Kcy * py * py
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
    std::vector<RealType> realRoots = poly.FindRealRoots();

    std::vector<RealType>::iterator ri;
    RealType r1, r2, alpha0;
    std::vector<std::pair<RealType,RealType> > rps;
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
        if (fabs(c1 + r1 * (b01 + r1 * a001) + r2 * (b11 + r2 * a111)) < 1e-6)
          { rps.push_back(std::make_pair(r1, r2)); }
        if (r1 > 1e-6) { //r1 non-negative
          r1 = -r1;
          if (fabs(c1 + r1 * (b01 + r1 * a001) + r2 * (b11 + r2 * a111)) <1e-6)
            { rps.push_back(std::make_pair(r1, r2)); }
        }
      }
    }
    // Consider combininig together the solving pair part w/ the searching
    // best solution part so that we don't need the pairs vector
    if (!rps.empty()) {
      RealType smallestDiff = HONKING_LARGE_VALUE;
      RealType diff;
      std::pair<RealType,RealType> bestPair = std::make_pair(1.0, 1.0);
      std::vector<std::pair<RealType,RealType> >::iterator rpi;
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
        std::cerr << "we choose r1 = " << bestPair.first
                  << " and r2 = " << bestPair.second << "\n";
#ifdef IS_MPI
      }
#endif

      RealType x, y, z;
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
      std::vector<StuntDouble*>::iterator sdi;
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
      exchangeSum_ += targetFlux_;
      //we may want to check whether the exchange has been successful
    } else {
      std::cerr << "exchange NOT performed!\n";//MPI incompatible
      failTrialCount_++;
    }

  }

  void RNEMD::doRNEMD() {

    switch(rnemdType_) {
    case rnemdKineticScale :
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

      if (rnemdLogWidth_ == midBin_ + 1)
        if (binNo > midBin_)
          binNo = nBins_ - binNo;

      RealType mass = sd->getMass();
      Vector3d vel = sd->getVel();
      RealType value;
      RealType xVal, yVal, zVal;

      switch(rnemdType_) {
      case rnemdKineticSwap :
      case rnemdKineticScale :
	
	value = mass * (vel[0]*vel[0] + vel[1]*vel[1] + 
			vel[2]*vel[2]);
	
        valueCount_[binNo] += 3;
	if (sd->isDirectional()) {
	  Vector3d angMom = sd->getJ();
	  Mat3x3d I = sd->getI();
          
	  if (sd->isLinear()) {
	    int i = sd->linearAxis();
	    int j = (i + 1) % 3;
	    int k = (i + 2) % 3;
	    value += angMom[j] * angMom[j] / I(j, j) + 
	      angMom[k] * angMom[k] / I(k, k);

            valueCount_[binNo] +=2;

	  } else {
	    value += angMom[0]*angMom[0]/I(0, 0) 
	      + angMom[1]*angMom[1]/I(1, 1) 
	      + angMom[2]*angMom[2]/I(2, 2);
            valueCount_[binNo] +=3;
	  }
	}
	value = value / PhysicalConstants::energyConvert / PhysicalConstants::kb;

	break;
      case rnemdPx :
      case rnemdPxScale :
	value = mass * vel[0];
        valueCount_[binNo]++;
        xVal = mass * vel.x() * vel.x() / PhysicalConstants::energyConvert
          / PhysicalConstants::kb;
        yVal = mass * vel.y() * vel.y() / PhysicalConstants::energyConvert
          / PhysicalConstants::kb;
        zVal = mass * vel.z() * vel.z() / PhysicalConstants::energyConvert
          / PhysicalConstants::kb;
        xTempHist_[binNo] += xVal;
        yTempHist_[binNo] += yVal;
        zTempHist_[binNo] += zVal;
	break;
      case rnemdPy :
      case rnemdPyScale :
	value = mass * vel[1];
        valueCount_[binNo]++;
	break;
      case rnemdPz :
      case rnemdPzScale :
	value = mass * vel[2];
        valueCount_[binNo]++;
	break;
      case rnemdUnknown : 
      default :
	break;
      }
      valueHist_[binNo] += value;
    }

  }

  void RNEMD::getStarted() {
    Snapshot* currentSnap_ = info_->getSnapshotManager()->getCurrentSnapshot();
    Stats& stat = currentSnap_->statData;
    stat[Stats::RNEMD_EXCHANGE_TOTAL] = exchangeSum_;
  }

  void RNEMD::getStatus() {

    Snapshot* currentSnap_ = info_->getSnapshotManager()->getCurrentSnapshot();
    Stats& stat = currentSnap_->statData;
    RealType time = currentSnap_->getTime();

    stat[Stats::RNEMD_EXCHANGE_TOTAL] = exchangeSum_;
    //or to be more meaningful, define another item as exchangeSum_ / time
    int j;

#ifdef IS_MPI

    // all processors have the same number of bins, and STL vectors pack their 
    // arrays, so in theory, this should be safe:

    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &valueHist_[0],
                              rnemdLogWidth_, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &valueCount_[0],
                              rnemdLogWidth_, MPI::INT, MPI::SUM);
    if (rnemdType_ == rnemdPx || rnemdType_ == rnemdPxScale) {
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &xTempHist_[0],
                                rnemdLogWidth_, MPI::REALTYPE, MPI::SUM);
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &yTempHist_[0],
                                rnemdLogWidth_, MPI::REALTYPE, MPI::SUM);
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &zTempHist_[0],
                                rnemdLogWidth_, MPI::REALTYPE, MPI::SUM);
    }
    // If we're the root node, should we print out the results
    int worldRank = MPI::COMM_WORLD.Get_rank();
    if (worldRank == 0) {
#endif
      rnemdLog_ << time;
      for (j = 0; j < rnemdLogWidth_; j++) {
        rnemdLog_ << "\t" << valueHist_[j] / (RealType)valueCount_[j];
      }
      rnemdLog_ << "\n";
      if (rnemdType_ == rnemdPx || rnemdType_ == rnemdPxScale ) {
        xTempLog_ << time;      
        for (j = 0; j < rnemdLogWidth_; j++) {
          xTempLog_ << "\t" << xTempHist_[j] / (RealType)valueCount_[j];
        }
        xTempLog_ << "\n";
        yTempLog_ << time;
        for (j = 0; j < rnemdLogWidth_; j++) {
          yTempLog_ << "\t" << yTempHist_[j] / (RealType)valueCount_[j];
        }
        yTempLog_ << "\n";
        zTempLog_ << time;
        for (j = 0; j < rnemdLogWidth_; j++) {
          zTempLog_ << "\t" << zTempHist_[j] / (RealType)valueCount_[j];
        }
        zTempLog_ << "\n";
      }
#ifdef IS_MPI
    }
#endif
    for (j = 0; j < rnemdLogWidth_; j++) {
      valueCount_[j] = 0;
      valueHist_[j] = 0.0;
    }
    if (rnemdType_ == rnemdPx || rnemdType_ == rnemdPxScale)
      for (j = 0; j < rnemdLogWidth_; j++) {
        xTempHist_[j] = 0.0;
        yTempHist_[j] = 0.0;
        zTempHist_[j] = 0.0;
      }
  }
}
