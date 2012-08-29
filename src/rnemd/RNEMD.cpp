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
#ifdef IS_MPI
#include <mpi.h>
#endif

#ifdef _MSC_VER
#define isnan(x) _isnan((x))
#define isinf(x) (!_finite(x) && !_isnan(x))
#endif

#define HONKING_LARGE_VALUE 1.0e10

using namespace std;
namespace OpenMD {
  
  RNEMD::RNEMD(SimInfo* info) : info_(info), evaluator_(info), seleMan_(info), 
                                usePeriodicBoundaryConditions_(info->getSimParams()->getUsePeriodicBoundaryConditions()) {

    trialCount_ = 0;
    failTrialCount_ = 0;
    failRootCount_ = 0;

    Globals * simParams = info->getSimParams();
    RNEMDParameters* rnemdParams = simParams->getRNEMDParameters();

    doRNEMD_ = rnemdParams->getUseRNEMD();
    if (!doRNEMD_) return;

    stringToMethod_["Swap"]  = rnemdSwap;
    stringToMethod_["NIVS"]  = rnemdNIVS;
    stringToMethod_["VSS"]   = rnemdVSS;

    stringToFluxType_["KE"]  = rnemdKE;
    stringToFluxType_["Px"]  = rnemdPx;
    stringToFluxType_["Py"]  = rnemdPy;
    stringToFluxType_["Pz"]  = rnemdPz;
    stringToFluxType_["Pvector"]  = rnemdPvector;
    stringToFluxType_["KE+Px"]  = rnemdKePx;
    stringToFluxType_["KE+Py"]  = rnemdKePy;
    stringToFluxType_["KE+Pvector"]  = rnemdKePvector;

    runTime_ = simParams->getRunTime();
    statusTime_ = simParams->getStatusTime();

    rnemdObjectSelection_ = rnemdParams->getObjectSelection();
    evaluator_.loadScriptString(rnemdObjectSelection_);
    seleMan_.setSelectionSet(evaluator_.evaluate());

    const string methStr = rnemdParams->getMethod();
    bool hasFluxType = rnemdParams->haveFluxType();

    string fluxStr;
    if (hasFluxType) {
      fluxStr = rnemdParams->getFluxType();
    } else {
      sprintf(painCave.errMsg, 
              "RNEMD: No fluxType was set in the md file.  This parameter,\n"
              "\twhich must be one of the following values:\n"
              "\tKE, Px, Py, Pz, Pvector, KE+Px, KE+Py, KE+Pvector\n"
              "\tmust be set to use RNEMD\n");
      painCave.isFatal = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    bool hasKineticFlux = rnemdParams->haveKineticFlux();
    bool hasMomentumFlux = rnemdParams->haveMomentumFlux();
    bool hasMomentumFluxVector = rnemdParams->haveMomentumFluxVector();
    bool hasSlabWidth = rnemdParams->haveSlabWidth();
    bool hasSlabACenter = rnemdParams->haveSlabACenter();
    bool hasSlabBCenter = rnemdParams->haveSlabBCenter();
    bool hasOutputFileName = rnemdParams->haveOutputFileName();
    bool hasOutputFields = rnemdParams->haveOutputFields();
    
    map<string, RNEMDMethod>::iterator i;
    i = stringToMethod_.find(methStr);
    if (i != stringToMethod_.end()) 
      rnemdMethod_ = i->second;
    else {
      sprintf(painCave.errMsg, 
              "RNEMD: The current method,\n"
              "\t\t%s is not one of the recognized\n"
              "\texchange methods: Swap, NIVS, or VSS\n",
              methStr.c_str());
      painCave.isFatal = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    map<string, RNEMDFluxType>::iterator j;
    j = stringToFluxType_.find(fluxStr);
    if (j != stringToFluxType_.end()) 
      rnemdFluxType_ = j->second;
    else {
      sprintf(painCave.errMsg, 
              "RNEMD: The current fluxType,\n"
              "\t\t%s\n"
              "\tis not one of the recognized flux types.\n",
              fluxStr.c_str());
      painCave.isFatal = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    bool methodFluxMismatch = false;
    bool hasCorrectFlux = false;
    switch(rnemdMethod_) {
    case rnemdSwap:
      switch (rnemdFluxType_) {
      case rnemdKE:
        hasCorrectFlux = hasKineticFlux;
        break;
      case rnemdPx:
      case rnemdPy:
      case rnemdPz:
        hasCorrectFlux = hasMomentumFlux;
        break;
      default :
        methodFluxMismatch = true;
        break;
      }
      break;
    case rnemdNIVS:
      switch (rnemdFluxType_) {
      case rnemdKE:
      case rnemdRotKE:
      case rnemdFullKE:
        hasCorrectFlux = hasKineticFlux;
        break;
      case rnemdPx:
      case rnemdPy:
      case rnemdPz:
        hasCorrectFlux = hasMomentumFlux;
        break;
      case rnemdKePx:
      case rnemdKePy:
        hasCorrectFlux = hasMomentumFlux && hasKineticFlux;
        break;
      default:
        methodFluxMismatch = true;
        break;
      }
      break;
    case rnemdVSS:
      switch (rnemdFluxType_) {
      case rnemdKE:
      case rnemdRotKE:
      case rnemdFullKE:
        hasCorrectFlux = hasKineticFlux;
        break;
      case rnemdPx:
      case rnemdPy:
      case rnemdPz:
        hasCorrectFlux = hasMomentumFlux;
        break;
      case rnemdPvector:
        hasCorrectFlux = hasMomentumFluxVector;
        break;
      case rnemdKePx:
      case rnemdKePy:
        hasCorrectFlux = hasMomentumFlux && hasKineticFlux;
        break;
      case rnemdKePvector:
        hasCorrectFlux = hasMomentumFluxVector && hasKineticFlux;
        break;
      default:
        methodFluxMismatch = true;
        break;
      }
    default:
      break;
    }

    if (methodFluxMismatch) {
      sprintf(painCave.errMsg, 
              "RNEMD: The current method,\n"
              "\t\t%s\n"
              "\tcannot be used with the current flux type, %s\n",
              methStr.c_str(), fluxStr.c_str());
      painCave.isFatal = 1;
      painCave.severity = OPENMD_ERROR;
      simError();        
    } 
    if (!hasCorrectFlux) {
      sprintf(painCave.errMsg, 
              "RNEMD: The current method, %s, and flux type, %s,\n"
              "\tdid not have the correct flux value specified. Options\n"
              "\tinclude: kineticFlux, momentumFlux, and momentumFluxVector\n",
              methStr.c_str(), fluxStr.c_str());
      painCave.isFatal = 1;
      painCave.severity = OPENMD_ERROR;
      simError();        
    } 

    if (hasKineticFlux) {
      // convert the kcal / mol / Angstroms^2 / fs values in the md file 
      // into  amu / fs^3:
      kineticFlux_ = rnemdParams->getKineticFlux() 
        * PhysicalConstants::energyConvert;
    } else {
      kineticFlux_ = 0.0;
    }
    if (hasMomentumFluxVector) {
      momentumFluxVector_ = rnemdParams->getMomentumFluxVector();
    } else {
      momentumFluxVector_ = V3Zero;
      if (hasMomentumFlux) {
        RealType momentumFlux = rnemdParams->getMomentumFlux();
        switch (rnemdFluxType_) {
        case rnemdPx:
          momentumFluxVector_.x() = momentumFlux;
          break;
        case rnemdPy:
          momentumFluxVector_.y() = momentumFlux;
          break;
        case rnemdPz:
          momentumFluxVector_.z() = momentumFlux;
          break;
        case rnemdKePx:
          momentumFluxVector_.x() = momentumFlux;
          break;
        case rnemdKePy:
          momentumFluxVector_.y() = momentumFlux;
          break;
        default:
          break;
        }
      }    
    }

    // do some sanity checking

    int selectionCount = seleMan_.getSelectionCount();
    int nIntegrable = info->getNGlobalIntegrableObjects();

    if (selectionCount > nIntegrable) {
      sprintf(painCave.errMsg, 
              "RNEMD: The current objectSelection,\n"
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

    areaAccumulator_ = new Accumulator();

    nBins_ = rnemdParams->getOutputBins();

    data_.resize(RNEMD::ENDINDEX);
    OutputData z;
    z.units =  "Angstroms";
    z.title =  "Z";
    z.dataType = "RealType";
    z.accumulator.reserve(nBins_);
    for (int i = 0; i < nBins_; i++) 
      z.accumulator.push_back( new Accumulator() );
    data_[Z] = z;
    outputMap_["Z"] =  Z;

    OutputData temperature;
    temperature.units =  "K";
    temperature.title =  "Temperature";
    temperature.dataType = "RealType";
    temperature.accumulator.reserve(nBins_);
    for (int i = 0; i < nBins_; i++) 
      temperature.accumulator.push_back( new Accumulator() );
    data_[TEMPERATURE] = temperature;
    outputMap_["TEMPERATURE"] =  TEMPERATURE;

    OutputData velocity;
    velocity.units = "angstroms/fs";
    velocity.title =  "Velocity";  
    velocity.dataType = "Vector3d";
    velocity.accumulator.reserve(nBins_);
    for (int i = 0; i < nBins_; i++) 
      velocity.accumulator.push_back( new VectorAccumulator() );
    data_[VELOCITY] = velocity;
    outputMap_["VELOCITY"] = VELOCITY;

    OutputData density;
    density.units =  "g cm^-3";
    density.title =  "Density";
    density.dataType = "RealType";
    density.accumulator.reserve(nBins_);
    for (int i = 0; i < nBins_; i++) 
      density.accumulator.push_back( new Accumulator() );
    data_[DENSITY] = density;
    outputMap_["DENSITY"] =  DENSITY;

    if (hasOutputFields) {
      parseOutputFileFormat(rnemdParams->getOutputFields());
    } else {
      outputMask_.set(Z);
      switch (rnemdFluxType_) {
      case rnemdKE:
      case rnemdRotKE:
      case rnemdFullKE:
        outputMask_.set(TEMPERATURE);
        break;
      case rnemdPx:
      case rnemdPy:
        outputMask_.set(VELOCITY);
        break;
      case rnemdPz:        
      case rnemdPvector:
        outputMask_.set(VELOCITY);
        outputMask_.set(DENSITY);
        break;
      case rnemdKePx:
      case rnemdKePy:
        outputMask_.set(TEMPERATURE);
        outputMask_.set(VELOCITY);
        break;
      case rnemdKePvector:
        outputMask_.set(TEMPERATURE);
        outputMask_.set(VELOCITY);
        outputMask_.set(DENSITY);        
        break;
      default:
        break;
      }
    }
      
    if (hasOutputFileName) {
      rnemdFileName_ = rnemdParams->getOutputFileName();
    } else {
      rnemdFileName_ = getPrefix(info->getFinalConfigFileName()) + ".rnemd";
    }          

    exchangeTime_ = rnemdParams->getExchangeTime();

    Snapshot* currentSnap_ = info->getSnapshotManager()->getCurrentSnapshot();
    Mat3x3d hmat = currentSnap_->getHmat();
  
    // Target exchange quantities (in each exchange) =  2 Lx Ly dt flux
    // Lx, Ly = box dimensions in x & y
    // dt = exchange time interval
    // flux = target flux

    RealType area = currentSnap_->getXYarea();
    kineticTarget_ = 2.0 * kineticFlux_ * exchangeTime_ * area;
    momentumTarget_ = 2.0 * momentumFluxVector_ * exchangeTime_ * area;

    // total exchange sums are zeroed out at the beginning:

    kineticExchange_ = 0.0;
    momentumExchange_ = V3Zero;

    if (hasSlabWidth) 
      slabWidth_ = rnemdParams->getSlabWidth();
    else
      slabWidth_ = hmat(2,2) / 10.0;
  
    if (hasSlabACenter) 
      slabACenter_ = rnemdParams->getSlabACenter();
    else 
      slabACenter_ = 0.0;
    
    if (hasSlabBCenter) 
      slabBCenter_ = rnemdParams->getSlabBCenter();
    else 
      slabBCenter_ = hmat(2,2) / 2.0;
    
  }
   
  RNEMD::~RNEMD() {
    if (!doRNEMD_) return;
#ifdef IS_MPI
    if (worldRank == 0) {
#endif

      writeOutputFile();

      rnemdFile_.close();
      
#ifdef IS_MPI
    }
#endif
  }
  
  bool RNEMD::inSlabA(Vector3d pos) {
    return (abs(pos.z() - slabACenter_) < 0.5*slabWidth_);
  }
  bool RNEMD::inSlabB(Vector3d pos) {
    return (abs(pos.z() - slabBCenter_) < 0.5*slabWidth_);
  }

  void RNEMD::doSwap() {
    if (!doRNEMD_) return;
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
      bool inA = inSlabA(pos);
      bool inB = inSlabB(pos);

      if (inA || inB) {
        
        RealType mass = sd->getMass();
        Vector3d vel = sd->getVel();
        RealType value;
        
        switch(rnemdFluxType_) {
        case rnemdKE :
          
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
        
        if (inA == 0) {
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
          
          switch(rnemdFluxType_) {
          case rnemdKE :
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
          
          switch(rnemdFluxType_) {
          case rnemdKE :
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
          
          switch(rnemdFluxType_) {
          case rnemdKE :
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
        
        switch(rnemdFluxType_) {
        case rnemdKE:
          kineticExchange_ += max_val - min_val;
          break;
        case rnemdPx:
          momentumExchange_.x() += max_val - min_val;
          break;
        case rnemdPy:
          momentumExchange_.y() += max_val - min_val;
          break;
        case rnemdPz:
          momentumExchange_.z() += max_val - min_val;
          break;
        default:
          break;
        }
      } else {        
        sprintf(painCave.errMsg, 
                "RNEMD::doSwap exchange NOT performed because min_val > max_val\n");
        painCave.isFatal = 0;
        painCave.severity = OPENMD_INFO;
        simError();        
        failTrialCount_++;
      }
    } else {
      sprintf(painCave.errMsg, 
              "RNEMD::doSwap exchange NOT performed because selected object\n"
	      "\twas not present in at least one of the two slabs.\n");
      painCave.isFatal = 0;
      painCave.severity = OPENMD_INFO;
      simError();        
      failTrialCount_++;
    }    
  }
  
  void RNEMD::doNIVS() {
    if (!doRNEMD_) return;
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
      bool inA = inSlabA(pos);
      bool inB = inSlabB(pos);

      if (inA || inB) {
               
        RealType mass = sd->getMass();
        Vector3d vel = sd->getVel();
       
        if (inA) {
          hotBin.push_back(sd);
          Phx += mass * vel.x();
          Phy += mass * vel.y();
          Phz += mass * vel.z();
          Khx += mass * vel.x() * vel.x();
          Khy += mass * vel.y() * vel.y();
          Khz += mass * vel.z() * vel.z();
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
        } else { 
          coldBin.push_back(sd);
          Pcx += mass * vel.x();
          Pcy += mass * vel.y();
          Pcz += mass * vel.z();
          Kcx += mass * vel.x() * vel.x();
          Kcy += mass * vel.y() * vel.y();
          Kcz += mass * vel.z() * vel.z();
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
    if ((rnemdFluxType_ == rnemdFullKE) ||
	(rnemdFluxType_ == rnemdRotKE)) {
      //may need sanity check Khw & Kcw > 0

      if (rnemdFluxType_ == rnemdFullKE) {
	c = 1.0 - kineticTarget_ / (Kcx + Kcy + Kcz + Kcw);
      } else {
	c = 1.0 - kineticTarget_ / Kcw;
      }

      if ((c > 0.81) && (c < 1.21)) {//restrict scaling coefficients
	c = sqrt(c);
	//std::cerr << "cold slab scaling coefficient: " << c << endl;
	//now convert to hotBin coefficient
	RealType w = 0.0;
	if (rnemdFluxType_ ==  rnemdFullKE) {
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
	    w = 1.0 + (kineticTarget_ 
                       + Khx * (1.0 - x * x) + Khy * (1.0 - y * y)
		       + Khz * (1.0 - z * z)) / Khw;
	  }//no need to calculate w if x, y or z is out of range
	} else {
	  w = 1.0 + kineticTarget_ / Khw;
	}
	if ((w > 0.81) && (w < 1.21)) {//restrict scaling coefficients
	  //if w is in the right range, so should be x, y, z.
	  vector<StuntDouble*>::iterator sdi;
	  Vector3d vel;
	  for (sdi = coldBin.begin(); sdi != coldBin.end(); sdi++) {
	    if (rnemdFluxType_ == rnemdFullKE) {
	      vel = (*sdi)->getVel() * c;
	      (*sdi)->setVel(vel);
	    }
	    if ((*sdi)->isDirectional()) {
	      Vector3d angMom = (*sdi)->getJ() * c;
	      (*sdi)->setJ(angMom);
	    }
	  }
	  w = sqrt(w);
	  // std::cerr << "xh= " << x << "\tyh= " << y << "\tzh= " << z
	  //           << "\twh= " << w << endl;
	  for (sdi = hotBin.begin(); sdi != hotBin.end(); sdi++) {
	    if (rnemdFluxType_ == rnemdFullKE) {
	      vel = (*sdi)->getVel();
	      vel.x() *= x;
	      vel.y() *= y;
	      vel.z() *= z;
	      (*sdi)->setVel(vel);
	    }
	    if ((*sdi)->isDirectional()) {
	      Vector3d angMom = (*sdi)->getJ() * w;
	      (*sdi)->setJ(angMom);
	    }
	  }
	  successfulScale = true;
	  kineticExchange_ += kineticTarget_;
	}
      }
    } else {
      RealType a000, a110, c0, a001, a111, b01, b11, c1;
      switch(rnemdFluxType_) {
      case rnemdKE :
	/* used hotBin coeff's & only scale x & y dimensions
	   RealType px = Phx / Pcx;
	   RealType py = Phy / Pcy;
	   a110 = Khy;
	   c0 = - Khx - Khy - kineticTarget_;
	   a000 = Khx;
	   a111 = Kcy * py * py;
	   b11 = -2.0 * Kcy * py * (1.0 + py);
	   c1 = Kcy * py * (2.0 + py) + Kcx * px * ( 2.0 + px) + kineticTarget_;
	   b01 = -2.0 * Kcx * px * (1.0 + px);
	   a001 = Kcx * px * px;
	*/
	//scale all three dimensions, let c_x = c_y
	a000 = Kcx + Kcy;
	a110 = Kcz;
	c0 = kineticTarget_ - Kcx - Kcy - Kcz;
	a001 = Khx * px * px + Khy * py * py;
	a111 = Khz * pz * pz;
	b01 = -2.0 * (Khx * px * (1.0 + px) + Khy * py * (1.0 + py));
	b11 = -2.0 * Khz * pz * (1.0 + pz);
	c1 = Khx * px * (2.0 + px) + Khy * py * (2.0 + py)
	  + Khz * pz * (2.0 + pz) - kineticTarget_;
	break;
      case rnemdPx :
	c = 1 - momentumTarget_.x() / Pcx;
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
      case rnemdPy :
	c = 1 - momentumTarget_.y() / Pcy;
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
      case rnemdPz ://we don't really do this, do we?
	c = 1 - momentumTarget_.z() / Pcz;
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
	  switch(rnemdFluxType_) {
	  case rnemdKE :
	    diff = fastpow(1.0 - r1, 2) + fastpow(1.0 - r2, 2)
	      + fastpow(r1 * r1 / r2 / r2 - Kcz/Kcx, 2)
	      + fastpow(r1 * r1 / r2 / r2 - Kcz/Kcy, 2);
	    break;
	  case rnemdPx :
	    diff = fastpow(1.0 - r1, 2) + fastpow(1.0 - r2, 2)
	      + fastpow(r1 * r1 / r2 / r2 - Kcz/Kcy, 2);
	    break;
	  case rnemdPy :
	    diff = fastpow(1.0 - r1, 2) + fastpow(1.0 - r2, 2)
	      + fastpow(r1 * r1 / r2 / r2 - Kcz/Kcx, 2);
	    break;
	  case rnemdPz :
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
	  // sprintf(painCave.errMsg, 
	  //         "RNEMD: roots r1= %lf\tr2 = %lf\n",
	  //         bestPair.first, bestPair.second);
	  // painCave.isFatal = 0;
	  // painCave.severity = OPENMD_INFO;
	  // simError();
#ifdef IS_MPI
	}
#endif
	
	switch(rnemdFluxType_) {
	case rnemdKE :
	  x = bestPair.first;
	  y = bestPair.first;
	  z = bestPair.second;
	  break;
	case rnemdPx :
	  x = c;
	  y = bestPair.first;
	  z = bestPair.second;
	  break;
	case rnemdPy :
	  x = bestPair.first;
	  y = c;
	  z = bestPair.second;
	  break;
	case rnemdPz :
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
        switch(rnemdFluxType_) {
        case rnemdKE :
          kineticExchange_ += kineticTarget_;
	  break;
	case rnemdPx :
	case rnemdPy :
	case rnemdPz :
          momentumExchange_ += momentumTarget_;
	  break;          
	default :
	  break;
	}       
      }
    }
    if (successfulScale != true) {
      sprintf(painCave.errMsg, 
              "RNEMD::doNIVS exchange NOT performed - roots that solve\n"
              "\tthe constraint equations may not exist or there may be\n"
              "\tno selected objects in one or both slabs.\n");
      painCave.isFatal = 0;
      painCave.severity = OPENMD_INFO;
      simError();        
      failTrialCount_++;
    }
  }

  void RNEMD::doVSS() {
    if (!doRNEMD_) return;
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
      bool inA = inSlabA(pos);
      bool inB = inSlabB(pos);
      
      if (inA || inB) {
        
        RealType mass = sd->getMass();
        Vector3d vel = sd->getVel();
       
        if (inA) {
          hotBin.push_back(sd);
	  //std::cerr << "before, velocity = " << vel << endl;
	  Ph += mass * vel;
	  //std::cerr << "after, velocity = " << vel << endl;
	  Mh += mass;
	  Kh += mass * vel.lengthSquare();
	  if (rnemdFluxType_ == rnemdFullKE) {
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
	  if (rnemdFluxType_ == rnemdFullKE) {
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
      Vector3d ac = -momentumTarget_ / Mc + vc;
      Vector3d acrec = -momentumTarget_ / Mc;
      RealType cNumerator = Kc - kineticTarget_ - 0.5 * Mc * ac.lengthSquare();
      if (cNumerator > 0.0) {
	RealType cDenominator = Kc - 0.5 * Mc * vc.lengthSquare();
	if (cDenominator > 0.0) {
	  RealType c = sqrt(cNumerator / cDenominator);
	  if ((c > 0.9) && (c < 1.1)) {//restrict scaling coefficients
	    Vector3d vh = Ph / Mh;
	    Vector3d ah = momentumTarget_ / Mh + vh;
	    Vector3d ahrec = momentumTarget_ / Mh;
	    RealType hNumerator = Kh + kineticTarget_
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
		    if (rnemdFluxType_ == rnemdFullKE) {
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
		    if (rnemdFluxType_ == rnemdFullKE) {
		      if ((*sdi)->isDirectional()) {
			Vector3d angMom = (*sdi)->getJ() * h;
			(*sdi)->setJ(angMom);
		      }
		    }
		  }
		  successfulExchange = true;
		  kineticExchange_ += kineticTarget_;
                  momentumExchange_ += momentumTarget_;
		}
	      }
	    }
	  }
	}
      }
    }
    if (successfulExchange != true) {
      sprintf(painCave.errMsg, 
              "RNEMD::doVSS exchange NOT performed - roots that solve\n"
              "\tthe constraint equations may not exist or there may be\n"
              "\tno selected objects in one or both slabs.\n");
      painCave.isFatal = 0;
      painCave.severity = OPENMD_INFO;
      simError();        
      failTrialCount_++;
    }
  }

  void RNEMD::doRNEMD() {
    if (!doRNEMD_) return;
    trialCount_++;
    switch(rnemdMethod_) {
    case rnemdSwap: 
      doSwap();
      break;
    case rnemdNIVS:
      doNIVS();
      break;
    case rnemdVSS:
      doVSS();
      break;
    case rnemdUnkownMethod:
    default :
      break;
    }
  }

  void RNEMD::collectData() {
    if (!doRNEMD_) return;
    Snapshot* currentSnap_ = info_->getSnapshotManager()->getCurrentSnapshot();
    Mat3x3d hmat = currentSnap_->getHmat();

    areaAccumulator_->add(currentSnap_->getXYarea());

    seleMan_.setSelectionSet(evaluator_.evaluate());

    int selei;
    StuntDouble* sd;
    int idx;

    vector<RealType> binMass(nBins_, 0.0);
    vector<RealType> binPx(nBins_, 0.0);
    vector<RealType> binPy(nBins_, 0.0);
    vector<RealType> binPz(nBins_, 0.0);
    vector<RealType> binKE(nBins_, 0.0);
    vector<int> binDOF(nBins_, 0);
    vector<int> binCount(nBins_, 0);

    // alternative approach, track all molecules instead of only those
    // selected for scaling/swapping:
    /*
    SimInfo::MoleculeIterator miter;
    vector<StuntDouble*>::iterator iiter;
    Molecule* mol;
    StuntDouble* sd;
    for (mol = info_->beginMolecule(miter); mol != NULL;
      mol = info_->nextMolecule(miter))
      sd is essentially sd
        for (sd = mol->beginIntegrableObject(iiter);
             sd != NULL;
             sd = mol->nextIntegrableObject(iiter))
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
      // Shift molecules by half a box to have bins start at 0
      // The modulo operator is used to wrap the case when we are 
      // beyond the end of the bins back to the beginning.
      int binNo = int(nBins_ * (pos.z() / hmat(2,2) + 0.5)) % nBins_;
     
      RealType mass = sd->getMass();
      Vector3d vel = sd->getVel();

      binCount[binNo]++;
      binMass[binNo] += mass;
      binPx[binNo] += mass*vel.x();
      binPy[binNo] += mass*vel.y();
      binPz[binNo] += mass*vel.z();
      binKE[binNo] += 0.5 * (mass * vel.lengthSquare());
      binDOF[binNo] += 3;

      if (sd->isDirectional()) {
        Vector3d angMom = sd->getJ();
        Mat3x3d I = sd->getI();
        if (sd->isLinear()) {
          int i = sd->linearAxis();
          int j = (i + 1) % 3;
          int k = (i + 2) % 3;
          binKE[binNo] += 0.5 * (angMom[j] * angMom[j] / I(j, j) + 
                                 angMom[k] * angMom[k] / I(k, k));
          binDOF[binNo] += 2;
        } else {
          binKE[binNo] += 0.5 * (angMom[0] * angMom[0] / I(0, 0) +
                                 angMom[1] * angMom[1] / I(1, 1) +
                                 angMom[2] * angMom[2] / I(2, 2));
          binDOF[binNo] += 3;
        }
      }
    }
    

#ifdef IS_MPI
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &binCount[0],
			      nBins_, MPI::INT, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &binMass[0],
			      nBins_, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &binPx[0],
			      nBins_, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &binPy[0],
			      nBins_, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &binPz[0],
			      nBins_, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &binKE[0],
			      nBins_, MPI::REALTYPE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &binDOF[0],
			      nBins_, MPI::INT, MPI::SUM);
#endif

    Vector3d vel;
    RealType den;
    RealType temp;
    RealType z;
    for (int i = 0; i < nBins_; i++) {
      z = (((RealType)i + 0.5) / (RealType)nBins_) * hmat(2,2);
      vel.x() = binPx[i] / binMass[i];
      vel.y() = binPy[i] / binMass[i];
      vel.z() = binPz[i] / binMass[i];

      den = binMass[i] * nBins_ * PhysicalConstants::densityConvert 
        / currentSnap_->getVolume() ;

      temp = 2.0 * binKE[i] / (binDOF[i] * PhysicalConstants::kb *
                               PhysicalConstants::energyConvert);

      for (unsigned int j = 0; j < outputMask_.size(); ++j) {
        if(outputMask_[j]) {
          switch(j) {
          case Z:
            (data_[j].accumulator[i])->add(z);
            break;
          case TEMPERATURE:
            data_[j].accumulator[i]->add(temp);
            break;
          case VELOCITY:
            dynamic_cast<VectorAccumulator *>(data_[j].accumulator[i])->add(vel);
            break;
          case DENSITY:
            data_[j].accumulator[i]->add(den);
            break;
          }
        }
      }
    }
  }

  void RNEMD::getStarted() {
    if (!doRNEMD_) return;
    collectData();
    writeOutputFile();
  }

  void RNEMD::parseOutputFileFormat(const std::string& format) {
    if (!doRNEMD_) return;
    StringTokenizer tokenizer(format, " ,;|\t\n\r");
    
    while(tokenizer.hasMoreTokens()) {
      std::string token(tokenizer.nextToken());
      toUpper(token);
      OutputMapType::iterator i = outputMap_.find(token);
      if (i != outputMap_.end()) {
        outputMask_.set(i->second);
      } else {
        sprintf( painCave.errMsg,
                 "RNEMD::parseOutputFileFormat: %s is not a recognized\n"
                 "\toutputFileFormat keyword.\n", token.c_str() );
        painCave.isFatal = 0;
        painCave.severity = OPENMD_ERROR;
        simError();            
      }
    }   
  }
  
  void RNEMD::writeOutputFile() {
    if (!doRNEMD_) return;
    
#ifdef IS_MPI
    // If we're the root node, should we print out the results
    int worldRank = MPI::COMM_WORLD.Get_rank();
    if (worldRank == 0) {
#endif
      rnemdFile_.open(rnemdFileName_.c_str(), std::ios::out | std::ios::trunc );
      
      if( !rnemdFile_ ){        
        sprintf( painCave.errMsg,
                 "Could not open \"%s\" for RNEMD output.\n",
                 rnemdFileName_.c_str());
        painCave.isFatal = 1;
        simError();
      }

      Snapshot* currentSnap_ = info_->getSnapshotManager()->getCurrentSnapshot();

      RealType time = currentSnap_->getTime();
      RealType avgArea;
      areaAccumulator_->getAverage(avgArea);
      RealType Jz = kineticExchange_ / (2.0 * time * avgArea) 
        / PhysicalConstants::energyConvert;
      Vector3d JzP = momentumExchange_ / (2.0 * time * avgArea);      

      rnemdFile_ << "#######################################################\n";
      rnemdFile_ << "# RNEMD {\n";

      map<string, RNEMDMethod>::iterator mi;
      for(mi = stringToMethod_.begin(); mi != stringToMethod_.end(); ++mi) {
        if ( (*mi).second == rnemdMethod_) 
          rnemdFile_ << "#    exchangeMethod  = \"" << (*mi).first << "\";\n";
      }
      map<string, RNEMDFluxType>::iterator fi;
      for(fi = stringToFluxType_.begin(); fi != stringToFluxType_.end(); ++fi) {
        if ( (*fi).second == rnemdFluxType_) 
          rnemdFile_ << "#    fluxType  = \"" << (*fi).first << "\";\n";
      }
      
      rnemdFile_ << "#    exchangeTime = " << exchangeTime_ << ";\n";

      rnemdFile_ << "#    objectSelection = \"" 
                 << rnemdObjectSelection_ << "\";\n";
      rnemdFile_ << "#    slabWidth = " << slabWidth_ << ";\n";
      rnemdFile_ << "#    slabAcenter = " << slabACenter_ << ";\n";
      rnemdFile_ << "#    slabBcenter = " << slabBCenter_ << ";\n";
      rnemdFile_ << "# }\n";
      rnemdFile_ << "#######################################################\n";
      rnemdFile_ << "# RNEMD report:\n";      
      rnemdFile_ << "#     running time = " << time << " fs\n";
      rnemdFile_ << "#     target flux:\n";
      rnemdFile_ << "#         kinetic = " 
                 << kineticFlux_ / PhysicalConstants::energyConvert 
                 << " (kcal/mol/A^2/fs)\n";
      rnemdFile_ << "#         momentum = " << momentumFluxVector_ 
                 << " (amu/A/fs^2)\n";
      rnemdFile_ << "#     target one-time exchanges:\n";
      rnemdFile_ << "#         kinetic = " 
                 << kineticTarget_ / PhysicalConstants::energyConvert 
                 << " (kcal/mol)\n";
      rnemdFile_ << "#         momentum = " << momentumTarget_ 
                 << " (amu*A/fs)\n";
      rnemdFile_ << "#     actual exchange totals:\n";
      rnemdFile_ << "#         kinetic = " 
                 << kineticExchange_ / PhysicalConstants::energyConvert 
                 << " (kcal/mol)\n";
      rnemdFile_ << "#         momentum = " << momentumExchange_ 
                 << " (amu*A/fs)\n";      
      rnemdFile_ << "#     actual flux:\n";
      rnemdFile_ << "#         kinetic = " << Jz
                 << " (kcal/mol/A^2/fs)\n";
      rnemdFile_ << "#         momentum = " << JzP 
                 << " (amu/A/fs^2)\n";
      rnemdFile_ << "#     exchange statistics:\n";
      rnemdFile_ << "#         attempted = " << trialCount_ << "\n";
      rnemdFile_ << "#         failed = " << failTrialCount_ << "\n";     
      if (rnemdMethod_ == rnemdNIVS) {
        rnemdFile_ << "#         NIVS root-check errors = "
                   << failRootCount_ << "\n";
      }
      rnemdFile_ << "#######################################################\n";
      
      
      
      //write title
      rnemdFile_ << "#";
      for (unsigned int i = 0; i < outputMask_.size(); ++i) {
        if (outputMask_[i]) {
          rnemdFile_ << "\t" << data_[i].title << 
            "(" << data_[i].units << ")";
          // add some extra tabs for column alignment
          if (data_[i].dataType == "Vector3d") rnemdFile_ << "\t\t";
        }
      }
      rnemdFile_ << std::endl;
      
      rnemdFile_.precision(8);
      
      for (int j = 0; j < nBins_; j++) {        
        
        for (unsigned int i = 0; i < outputMask_.size(); ++i) {
          if (outputMask_[i]) {
            if (data_[i].dataType == "RealType")
              writeReal(i,j);
            else if (data_[i].dataType == "Vector3d")
              writeVector(i,j);
            else {
              sprintf( painCave.errMsg,
                       "RNEMD found an unknown data type for: %s ",
                       data_[i].title.c_str());
              painCave.isFatal = 1;
              simError();
            }
          }
        }
        rnemdFile_ << std::endl;
        
      }        

      rnemdFile_ << "#######################################################\n";
      rnemdFile_ << "# Standard Deviations in those quantities follow:\n";
      rnemdFile_ << "#######################################################\n";


      for (int j = 0; j < nBins_; j++) {        
        rnemdFile_ << "#";
        for (unsigned int i = 0; i < outputMask_.size(); ++i) {
          if (outputMask_[i]) {
            if (data_[i].dataType == "RealType")
              writeRealStdDev(i,j);
            else if (data_[i].dataType == "Vector3d")
              writeVectorStdDev(i,j);
            else {
              sprintf( painCave.errMsg,
                       "RNEMD found an unknown data type for: %s ",
                       data_[i].title.c_str());
              painCave.isFatal = 1;
              simError();
            }
          }
        }
        rnemdFile_ << std::endl;
        
      }        
      
      rnemdFile_.flush();
      rnemdFile_.close();
      
#ifdef IS_MPI
    }
#endif
    
  }
  
  void RNEMD::writeReal(int index, unsigned int bin) {
    if (!doRNEMD_) return;
    assert(index >=0 && index < ENDINDEX);
    assert(bin < nBins_);
    RealType s;
    
    data_[index].accumulator[bin]->getAverage(s);
    
    if (! isinf(s) && ! isnan(s)) {
      rnemdFile_ << "\t" << s;
    } else{
      sprintf( painCave.errMsg,
               "RNEMD detected a numerical error writing: %s for bin %d",
               data_[index].title.c_str(), bin);
      painCave.isFatal = 1;
      simError();
    }    
  }
  
  void RNEMD::writeVector(int index, unsigned int bin) {
    if (!doRNEMD_) return;
    assert(index >=0 && index < ENDINDEX);
    assert(bin < nBins_);
    Vector3d s;
    dynamic_cast<VectorAccumulator*>(data_[index].accumulator[bin])->getAverage(s);
    if (isinf(s[0]) || isnan(s[0]) || 
        isinf(s[1]) || isnan(s[1]) || 
        isinf(s[2]) || isnan(s[2]) ) {      
      sprintf( painCave.errMsg,
               "RNEMD detected a numerical error writing: %s for bin %d",
               data_[index].title.c_str(), bin);
      painCave.isFatal = 1;
      simError();
    } else {
      rnemdFile_ << "\t" << s[0] << "\t" << s[1] << "\t" << s[2];
    }
  }  

  void RNEMD::writeRealStdDev(int index, unsigned int bin) {
    if (!doRNEMD_) return;
    assert(index >=0 && index < ENDINDEX);
    assert(bin < nBins_);
    RealType s;
    
    data_[index].accumulator[bin]->getStdDev(s);
    
    if (! isinf(s) && ! isnan(s)) {
      rnemdFile_ << "\t" << s;
    } else{
      sprintf( painCave.errMsg,
               "RNEMD detected a numerical error writing: %s std. dev. for bin %d",
               data_[index].title.c_str(), bin);
      painCave.isFatal = 1;
      simError();
    }    
  }
  
  void RNEMD::writeVectorStdDev(int index, unsigned int bin) {
    if (!doRNEMD_) return;
    assert(index >=0 && index < ENDINDEX);
    assert(bin < nBins_);
    Vector3d s;
    dynamic_cast<VectorAccumulator*>(data_[index].accumulator[bin])->getStdDev(s);
    if (isinf(s[0]) || isnan(s[0]) || 
        isinf(s[1]) || isnan(s[1]) || 
        isinf(s[2]) || isnan(s[2]) ) {      
      sprintf( painCave.errMsg,
               "RNEMD detected a numerical error writing: %s std. dev. for bin %d",
               data_[index].title.c_str(), bin);
      painCave.isFatal = 1;
      simError();
    } else {
      rnemdFile_ << "\t" << s[0] << "\t" << s[1] << "\t" << s[2];
    }
  }  
}

