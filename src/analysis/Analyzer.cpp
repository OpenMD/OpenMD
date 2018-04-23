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
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */
#ifdef IS_MPI
#include <mpi.h>
#endif

#include <cmath>
#include <sstream>
#include <string>

#include "analysis/Analyzer.hpp"
#include "math/Vector3.hpp"
#include "math/Vector.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Polynomial.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/StuntDouble.hpp"
#include "utils/Constants.hpp"
#include "utils/Tuple.hpp"
#include "brains/Thermo.hpp"
#include "math/ConvexHull.hpp"

#include "analysis/StaticAnalyser.hpp"
#include "analysis/GofR.hpp"
#include "analysis/GofZ.hpp"
#include "analysis/GofRZ.hpp"
#include "analysis/GofRAngle.hpp"
#include "analysis/GofAngle2.hpp"
#include "analysis/GofRAngle2.hpp"
#include "analysis/GofXyz.hpp"
#include "analysis/TwoDGofR.hpp"
#include "analysis/P2OrderParameter.hpp"
#include "analysis/BondOrderParameter.hpp"
#include "analysis/BOPofR.hpp"
#include "analysis/RippleOP.hpp"
#include "analysis/SCDOrderParameter.hpp"
#include "analysis/DensityPlot.hpp"
#include "analysis/ObjectCount.hpp"
#include "analysis/RhoZ.hpp"
#include "analysis/PipeDensity.hpp"
#include "analysis/pAngle.hpp"
#include "analysis/BondAngleDistribution.hpp"
#include "analysis/NanoVolume.hpp"
#include "analysis/NanoLength.hpp"
#if defined(HAVE_FFTW_H) || defined(HAVE_DFFTW_H) || defined(HAVE_FFTW3_H)
#include "analysis/Hxy.hpp"
#endif
#include "analysis/RhoR.hpp"
#include "analysis/AngleR.hpp"
#include "analysis/TetrahedralityParam.hpp"
#include "analysis/TetrahedralityParamZ.hpp"
#include "analysis/TetrahedralityParamXYZ.hpp"
#include "analysis/TetrahedralityParamDens.hpp"
#include "analysis/RNEMDStats.hpp"
#include "analysis/NitrileFrequencyMap.hpp"
#include "analysis/MultipoleSum.hpp"
#include "analysis/CoordinationNumber.hpp"
#include "analysis/HBondGeometric.hpp"
#include "analysis/PotDiff.hpp"
#include "analysis/TetrahedralityHBMatrix.hpp"
#include "analysis/Kirkwood.hpp"
#include "analysis/Field.hpp"
#include "analysis/VelocityZ.hpp"


#ifdef _MSC_VER
#define isnan(x) _isnan((x))
#define isinf(x) (!_finite(x) && !_isnan(x))
#else
#define isnan(x) std::isnan((x))
#define isinf(x) std::isinf((x))
#endif

#define HONKING_LARGE_VALUE 1.0e10

using namespace std;
namespace OpenMD {
  
  Analyzer::Analyzer(SimInfo* info) : info_(info),evaluator_(info_),
				      seleMan_(info_){

    Globals* simParams = info->getSimParams();
    AnalyzerParameters* analyzerParams = simParams->getAnalyzerParameters();
    
    doAnalyzer_ = analyzerParams->getUseAnalyzer();
    if (!doAnalyzer_) return;

    stringToMethod_["bo"] = analyzerBo;
    stringToMethod_["ior"] = analyzerIor;
    stringToMethod_["for"] = analyzerFor;
    stringToMethod_["bad"] = analyzerBad;
    stringToMethod_["count"] = analyzerCount;
    stringToMethod_["gofr"] = analyzerGofr;
    stringToMethod_["gofz"] = analyzerGofz;
    stringToMethod_["rTheta"] = analyzerRTheta;
    stringToMethod_["rOmega"] = analyzerROmega;
    stringToMethod_["rz"] = analyzerRz;
    stringToMethod_["thetaOmega"] = analyzerThetaOmega;
    stringToMethod_["rThetaOmega"] = analyzerRThetaOmega;
    stringToMethod_["gxyz"] = analyzerGxyz;
    stringToMethod_["twodgofr"] = analyzerTwoDGofr;
    stringToMethod_["p2"] = analyzerP2;
    stringToMethod_["scd"] = analyzerSCD;
    stringToMethod_["density"] = analyzerDensity;
    stringToMethod_["slabDensity"] = analyzerSlabDensity;
    stringToMethod_["pipeDensity"] = analyzerPipeDensity;
    stringToMethod_["pAngle"] = analyzerPAngle;
    stringToMethod_["hxy"] = analyzerHxy;
    stringToMethod_["rhoR"] = analyzerRhoR;
    stringToMethod_["angleR"] = analyzerAngleR;
    stringToMethod_["hullVol"] = analyzerHullVol;
    stringToMethod_["rodLength"] = analyzerRodLength;
    stringToMethod_["tetParam"] = analyzerTetParam;
    stringToMethod_["tetParamZ"] = analyzerTetParamZ;
    stringToMethod_["tetParamDens"] = analyzerTetParamDens;
    stringToMethod_["tetParamXYZ"] = analyzerTetParamXYZ;
    stringToMethod_["rnemdZ"] = analyzerRNEMDz;
    stringToMethod_["rnemdR"] = analyzerRNEMDr;
    stringToMethod_["rnemdRT"] = analyzerRNEMDrt;
    stringToMethod_["nitrile"] = analyzerNitrile;
    stringToMethod_["multipole"] = analyzerMultipole;
    stringToMethod_["surfDiffusion"] = analyzerSurfDiffusion;
    stringToMethod_["cn"] = analyzerCN;
    stringToMethod_["scn"] = analyzerSCN;
    stringToMethod_["gcn"] = analyzerGCN;
    stringToMethod_["hBond"] = analyzerHBond;
    stringToMethod_["potDiff"] = analyzerPotDiff;
    stringToMethod_["tetHB"] = analyzerTetHB;
    stringToMethod_["kirkwood"] = analyzerKirkwood;
    stringToMethod_["kirkwoodQ"] = analyzerKirkwoodQ;
    stringToMethod_["densityField"] = analyzerDensityField;
    stringToMethod_["velocityField"] = analyzerVelocityField;
    stringToMethod_["velocityZ"] = analyzerVelocityZ;
    
    const string methStr = analyzerParams->getMethod();

    map<string, AnalyzerMethod>::iterator i;
    i = stringToMethod_.find(methStr);
    if (i != stringToMethod_.end()) 
      analyzerMethod_ = i->second;
    else {
      sprintf(painCave.errMsg, 
              "Analyzer: The current method,\n"
              "\t\t%s is not one of the recognized\n"
              "\tanalysis methods:\n"
	      "\tbo, ior, for, bad, count, gofr, gofz, rTheta, rOmega\n"
	      "\trz, thetaOmega, rThetaOmega, gxyz, twodgofr, p2, scd\n"
	      "\tdensity, slabDensity, pipeDensity, pAngle, hxy, rhoR, angleR\n"
	      "\thullVol, rodLength, tetParam, tetParamZ, tetParamDens\n"
	      "\ttetParamXYZ, rnemdZ, rnemdR, rnemdRT, nitrile, multipole \n"
	      "\tsurfDiffusion, cn, scn, gcn, hBond, potDiff, tetHB \n"
	      "\tkirkwood, kirkwoodQ, densityField, velocityField, velocityZ\n",
              methStr.c_str());
      painCave.isFatal = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    // Test for all of the parameters now, this will prevent redundant code
    // when we call analyses.
    bool hasStep = analyzerParams->haveStep();
    bool hasNBins = analyzerParams->haveNBins();
    bool hasNBinsX = analyzerParams->haveNBinsX();
    bool hasNBinsY = analyzerParams->haveNBinsY();
    bool hasNBinsZ = analyzerParams->haveNBinsZ();
    bool hasNRBins = analyzerParams->haveNRBins();
    bool hasNAngleBins = analyzerParams->haveNAngleBins();
    bool hasRCut = analyzerParams->haveRCut();
    bool hasOOCut = analyzerParams->haveOOCut();
    bool hasThetaCut = analyzerParams->haveThetaCut();
    bool hasOHCut = analyzerParams->haveOHCut();
    bool hasDz = analyzerParams->haveDz();
    bool hasLength = analyzerParams->haveLength();
    bool hasZLength = analyzerParams->haveZLength();
    bool hasZOffSet = analyzerParams->haveZOffSet();
    bool hasSele1 = analyzerParams->haveSele1();
    bool hasSele2 = analyzerParams->haveSele2();
    bool hasSele3 = analyzerParams->haveSele3();
    bool hasComSele = analyzerParams->haveComSele();
    bool hasSeleOffSet = analyzerParams->haveSeleOffSet();
    bool hasSeleOffSet2 = analyzerParams->haveSeleOffSet2();
    bool hasMolName = analyzerParams->haveMolName();
    bool hasBegin = analyzerParams->haveBegin();
    bool hasEnd = analyzerParams->haveEnd();
    bool hasRadius = analyzerParams->haveRadius();
    bool hasVoxelSize = analyzerParams->haveVoxelSize();
    bool hasGaussWidth = analyzerParams->haveGaussWidth();
    bool hasPrivilegedAxis = analyzerParams->havePrivilegedAxis();
    bool hasPrivilegedAxis2 = analyzerParams->havePrivilegedAxis2();

    // Store parameters for analyses, again prevent redudant code.
    if (hasStep)
      step = analyzerParams->getStep();
    if (hasNBins)
      nbins = analyzerParams->getNBins();
    if (hasNBinsX)
      nbinsX = analyzerParams->getNBinsX();
    if (hasNBinsY)
      nbinsY = analyzerParams->getNBinsY();
    if (hasNBinsZ)
      nbinsZ = analyzerParams->getNBinsZ();
    if (hasNRBins)
      nrbins = analyzerParams->getNRBins();
    if (hasNAngleBins)
      nAngleBins = analyzerParams->getNAngleBins();
    if (hasRCut)
      rCut = analyzerParams->getRCut();
    if (hasOOCut)
      ooCut = analyzerParams->getOOCut();
    if (hasThetaCut)
      thetaCut = analyzerParams->getThetaCut();
    if (hasOHCut)
      OHCut = analyzerParams->getOHCut();
    if (hasDz)
      dz = analyzerParams->getDz();
    if (hasLength)
      length = analyzerParams->getLength();
    if (hasZLength)
      zLength = analyzerParams->getZLength();
    if (hasZOffSet)
      zOffSet = analyzerParams->getZOffSet();
    if (hasSele1)
      sele1 = analyzerParams->getSele1();
    if (hasSele2)
      sele2 = analyzerParams->getSele2();
    if (hasSele3)
      sele3 = analyzerParams->getSele3();
    if (hasComSele)
      comSele = analyzerParams->getComSele();
    if (hasSeleOffSet)
      seleOffSet = analyzerParams->getSeleOffSet();
    if (hasSeleOffSet2)
      seleOffSet2 = analyzerParams->getSeleOffSet2();
    if (hasMolName)
      molName = analyzerParams->getMolName();
    if (hasBegin)
      begin = analyzerParams->getBegin();
    if (hasEnd)
      end = analyzerParams->getEnd();
    if (hasRadius)
      radius = analyzerParams->getRadius();
    if (hasVoxelSize)
      voxelSize = analyzerParams->getVoxelSize();
    if (hasGaussWidth)
      gaussWidth = analyzerParams->getGaussWidth();
    if (hasPrivilegedAxis)
      privilegedAxis = analyzerParams->getPrivilegedAxis();
    if (hasPrivilegedAxis2)
      privilegedAxis2 = analyzerParams->getPrivilegedAxis2();

    
    // convert privilegedAxis to corresponding integer
    // x axis -> 0
    // y axis -> 1
    // z axis -> 2 (default)
    
    // if (hasPrivilegedAxis){
    //   std::string privilegedAxisStr = analyzerParams->getPrivilegedAxis();
    //   int privilegedAxis;
    //   switch (privilegedAxisStr) {
    //   case 'x':
    // 	privilegedAxis = 0;
    // 	break;
    //   case 'y':
    // 	privilegedAxis = 1;
    // 	break;
    //   case 'z':
    //   default:
    // 	privilegedAxis = 2;
    // 	break;
    //   }
    // }

    // if (hasPrivilegedAxis2){
    //   std::string privilegedAxis2Str = analyzerParams->getPrivilegedAxis2();
    //   int privilegedAxis2;
    //   switch (privilegedAxis2Str) {
    //   case 'x':
    // 	privilegedAxis2 = 0;
    // 	break;
    //   case 'y':
    // 	privilegedAxis2 = 1;
    // 	break;
    //   case 'z':
    //   default:
    // 	privilegedAxis2 = 2;
    // 	break;
    //   }
    // }

    if (!hasLength) {
      hmat = info->getSnapshotManager()->getCurrentSnapshot()->getHmat();
      maxLen = std::min(std::min(hmat(0, 0), hmat(1, 1)), hmat(2, 2)) /2.0;
      zmaxLen = hmat(2,2);
      hasLength = true;
    }
    
    
    switch(analyzerMethod_) {
    case analyzerBo:
      if (hasSele1 && hasRCut && hasNBins) {
	analyser_ = new BondOrderParameter(info, sele1, rCut, nbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'bo' requires the following specified\n"
		"\tSele1, Sele2, Length, and NRBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerIor:
      if (hasSele1 && hasSele2 && hasLength && hasNBins && hasPrivilegedAxis) {
	analyser_ = new IcosahedralOfR(info, sele1, rCut, nrbins, maxLen);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'ior' requires the following specified\n"
		"\tSele1, Sele2, Length, NBins, privilegedAxis\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerFor:
      if (hasSele1 && hasRCut && hasNRBins && hasLength && hasPrivilegedAxis) {
	analyser_ = new FCCOfR(info, sele1, rCut, nrbins, length);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'for' requires the follwing specified\n"
		"\tSele1, rCut, nRBins, Length\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerBad:
      std::cerr << "in analyzerBad" << endl;
      if (hasSele1 && hasRCut && hasNBins) {
	analyser_ = new BondAngleDistribution(info, sele1, rCut, nbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'bad' requires the follwing specified\n"
		"\tSele1, rCut, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerCount:
      std::cerr << "in analyzerCount" << endl;
      if (hasSele1) {
	analyser_ = new ObjectCount(info, sele1);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'count' requires the follwing specified\n"
		"\tSele1\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerGofr:
      if (hasSele1 && hasSele2 && hasLength && hasNRBins) {
	analyser_ = new GofR(info, sele1, sele2, length, nrbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'gofr' requires the follwing specified\n"
		"\tSele1, Sele2, Length, nRBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerGofz:
      if (hasSele1 && hasSele2 && hasLength && hasNBins && hasPrivilegedAxis) {
	analyser_ = new GofZ(info, sele1, sele2, length, nbins, privilegedAxis);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'gofz' requires the follwing specified\n"
		"\tSele1, Sele2, Length, NBins, PrivilegedAxis\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerRTheta:
      if (hasSele1 && hasSele2 && hasSele3 && hasLength && hasNRBins
	  && hasNAngleBins) {
	analyser_ = new GofRTheta(info, sele1, sele2, sele3, length, nrbins, nAngleBins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'rtheta' requires the follwing specified\n"
		"\tSele1, Sele2, Sele3, nAngleBins, nRBins, Length\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerROmega:
      if (hasSele1 && hasSele2 && hasSele3 && hasLength && hasNRBins
	  && hasNAngleBins) {
	analyser_ = new GofROmega(info, sele1, sele2, sele3, length, nrbins, nAngleBins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'romega' requires the follwing specified\n"
		"\tSele1, Sele2, Sele3, nAngleBins, nRBins, Length\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerRz:
      if (hasSele1 && hasSele2 && hasLength && hasNRBins && hasNBinsZ
	  && hasPrivilegedAxis) {
	analyser_ = new GofRZ(info, sele1, sele2, length, maxLength, nrbins, nbinsZ, privilegedAxis);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'rz' requires the follwing specified\n"
		"\tSele1, Sele2, Length, NRBins, NBinsZ, PrivilegedAxis\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerThetaOmega:
      if (hasSele1 && hasSele2 && hasSele3 && hasNAngleBins) {
	anlyser_ = new GofAngle2(info, sele1, sele2, sele3, nagnlebins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'thetaOmega' requires the follwing specified\n"
		"\tSele1, Sele2, Sele3, nAngleBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerRThetaOmega:
      if (hasSele1 && hasSele2 && hasSele3 && hasNRBins && hasNAngleBins) {
	analyser_ = new GofRAngle2(info, sele1, sele2, sele3, maxLength, nrbins, nAngleBins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'rThetaOmega' requires the follwing specified\n"
		"\tSele1, Sele2, Sele3, nRBins, nAngleBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerGxyz:
      if (hasSele1 && hasSele2 && hasNBins) {
	analyser_ = new GofXyz(info, sele1, sele2, maxLength, nbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'gxyz' requires the follwing specified\n"
		"\tSele1, Sele2, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerTwoDGofr:
      if (hasSele1 && hasSele2 && hasDz && hasNRBins) {
	analyser_ = new TwoDGofR(indo, sele1, sele2, maxLength, dz, nrbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'twodgofr' requires the follwing specified\n"
		"\tSele1, rCut, nRBins, Length\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerP2:
      if (hasSele1) {
	if (hasSele2) {
	  analyser_ = new P2OrderParameter(info, sele1, sele2);
	  break;
	}
	else {
	  if (hasSeleOffSet) {
	    analyser_ = new P2OrderParameter(info, sele1, seleOffSet);
	    break;
	  }
	  else {
	    analyser_ = new P2OrderParameter(info, sele1);
	    break;
	  }
	}
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'p2' requires the follwing specified\n"
		"\tSele1\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerSCD:
      if (hasSele1 && hasSele2 && hasSele3) {
	analyser_ = new SCDOrderParameter(info, sele1, sele2, sele3);
	break;
      } else {
		sprintf(painCave.errMsg,
		"Analyzer: 'scd' requires the follwing specified\n"
		"\tSele1, Sele2, Sele3\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerDensity:
      if (hasSele1 && hasSele2 && hasNBins) {
	analyser_ = new DensityPlot(info, sele1, sele2, maxLength, nbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'density' requires the follwing specified\n"
		"\tSele1, Sele2, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerSlabDensity:
      if (hasSele1 && hasNBins && hasPrivilegedAxis) {
	analyser_ = new RhoZ(info, sele1, nbins, privilegedAxis);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'slabDensity' requires the follwing specified\n"
		"\tSele1, Sele2, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerPipeDensity:
    case analyzerPAngle:
      if (hasSele1) {
	if (hasSele2) {
	  analyser_ = new pAngle(info, sele1, sele2, nbins);
	  break;
	} else {
	  if (hasSeleOffSet) {
	    if (hasSeleOffSet2) {
	      analyser_ = new pAngle(info, sele1, seleOffSet, seleOffSet2, nbins);
	      break;
	    } else {
	      analyser_ = new pAngle(info, sele1, seleOffSet, nbins);
	      break;
	    }
	  }
	  analyser_ = new pAngle(info, sele1, nbins);
	  break;
	}
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'pAngle' requires the follwing specified\n"
		"\tSele1, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerHxy:
#if defined(HAVE_FFTW_H) || defined(HAVE_DFFTW_H) || defined(HAVE_FFTW3_H)
      if (hasSele1 && hasNBinsX && hasNBinsY && hasNBinsZ && hasNBins) {
	analyser_ = new Hxy(info, sele1, nbinsX, nbinsY, nbinsZ, nbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'hxy' requires the follwing specified\n"
		"\tSele1, NBinsX, NBinsY, NBinsZ, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }	
#endif
    case analyzerRhoR:
      if (hasSele1 && hasNRBins && hasRadius) {
	analyser_ = new RhoR(info, sele1, maxLength, nrbins, radius);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'rhor' requires the follwing specified\n"
		"\tSele1, nRBins, Radius\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerAngleR:
      if (hasSele1 && hasNRBins) {
	analyser_ = new AngleR(info, sele1, maxLength, nrbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'angleR' requires the follwing specified\n"
		"\tSele1, nRBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerHullVol:
      if (hasSele1) {
	analyser_ = new NanoVolume(info, sele1);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'hullVol' requires the follwing specified\n"
		"\tSele1\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerRodLength:
      if (hasSele1) {
	analyser_ = new NanoLength(info, sele1);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'rodLength' requires the follwing specified\n"
		"\tSele1\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerTetParam:
      if (hasSele1 && hasRCut && hasNBins) {
	analyser_ = new TetrahedralityParam(info, sele1, rCut, nbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'tetParam' requires the follwing specified\n"
		"\tSele1, rCut, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerTetParamZ:
      if (hasSele1 && hasSele2 && hasRCut && hasNBins && hasPrivilegedAxis) {
	analyser_ = new TetrahedralityParamZ(info, sele1, sele2, rCut, nbins, privilegedAxis);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'tetParamZ' requires the follwing specified\n"
		"\tSele1, Sele2, rCut, nBins, PrivilegedAxis\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerTetParamDens:
      if (hasSele1 && hasSele2 && hasRCut && hasNBins) {
	analyser_ = new TetrahedralityParamDens(info, sele1, sele2, rCut, nbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'tetParamDens' requires the follwing specified\n"
		"\tSele1, Sele2, RCut, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerTetParamXYZ:
      if (hasSele1 && hasSele2 && hasRCut && hasVoxelSize && hasGaussWidth) {
	analyser_ = new TetrahedralityParamXYZ(info, sele1, sele2, rCut, voxelSize, gaussWidth);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'pAngle' requires the follwing specified\n"
		"\tSele1, Sele2, RCut, VoxelSize, GaussWidth, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerRNEMDz:
      if (hasSele1 && hasNBins && hasPrivilegedAxis) {
	analyser_ = new RNEMDZ(info, sele1, nbins, privilegedAxis);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'rnemdz' requires the follwing specified\n"
		"\tSele1, nBins, PrivilegedAxis\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerRNEMDr:
      if (hasSele1 && hasNRBins) {
	analyser_ = new RNEMDR(info, sele1, nrbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'rnemdR' requires the follwing specified\n"
		"\tSele1, nRBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerRNEMDrt:
      if (hasSele1 && hasNRBins && hasNAngleBins) {
	analyser_ = new RNEMDRTheta(info, sele1, nrbins, nAngleBins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'rnemdRT' requires the follwing specified\n"
		"\tSele1, nRBins, nAngleBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerNitrile:
      if (hasSele1 && hasNBins) {
	analyser_ = new NitrileFrequencyMap(info, sele1, nbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'nitrile' requires the follwing specified\n"
		"\tSele1, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerMultipole:
      if (hasSele1 && hasNBins) {
	analyser_ = new MultipoleSum(info, sele1, maxLength, nbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'multipole' requires the follwing specified\n"
		"\tSele1, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerSurfDiffusion:
      if (hasSele1) {
	analyser_ = new SurfaceDiffusion(info, sele1, maxLength);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'surfDiffusion' requires the follwing specified\n"
		"\tSele1\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerCN:
      if (hasSele1 && hasSele2 && hasRCut && hasNBins) {
	analyser_ = new CoordinationNumber(info, sele1, sele2, rCut, nbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'cn' requires the follwing specified\n"
		"\tSele1, Sele2, RCut, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerSCN:
      if (hasSele1 && hasSele2 && hasRCut && hasNBins) {
	analyser_ = new SCN(info, sele1, sele2, rCut, nbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'scn' requires the follwing specified\n"
		"\tSele1, Sele2, RCut, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerGCN:
      if (hasSele1 && hasSele2 && hasRCut && hasNBins) {
	analyser_ = new GCN(info, sele1, sele2, rCut, nbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'gcn' requires the follwing specified\n"
		"\tSele1, Sele2, RCut, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerHBond:
      if (hasSele1 && hasSele2 && hasRCut && hasThetaCut && hasNBins) {
	analyser_ = new HBondGeometric(info, sele1, sele2, rCut, thetaCut, nbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'hbond' requires the follwing specified\n"
		"\tSele1, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerPotDiff:
      if (hasSele1) {
	analyser_ = new PotDiff(info, sele1);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'potDiff' requires the follwing specified\n"
		"\tSele1, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerTetHB:
      if (hasSele1 && hasRCut && hasOOCut && hasThetaCut && hasOHCut && hasNBins) {
	analyser_ = new TetrahedralityHBMatrix(info, sele1, rCut, ooCut, thetaCut, OHCut, nbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'hbtet' requires the follwing specified\n"
		"\tSele1, RCut, OOCut, ThetaCut, OHCut, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerKirkwood:
      if (hasSele1 && hasSele2 && hasNRBins) {
	analyser_ = new Kirkwood(info, sele1, sele2, maxLength, nrbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'kirkwood' requires the follwing specified\n"
		"\tSele1, Sele2, nRBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerKirkwoodQ:
      if (hasSele1 && hasSele2 && hasNRBins) {
	analyser_ = new KirkwoodQ(info, sele1, sele2, maxLength, nrbins);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'kirkwoodQ' requires the follwing specified\n"
		"\tSele1, Sele2, nRBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerDensityField:
      if (hasSele1 && hasVoxelSize) {
	analyser_ = new DensityField(info, sele1, voxelSize);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'densityField' requires the follwing specified\n"
		"\tSele1, VoxelSize\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerVelocityField:
      if (hasSele1 && hasVoxelSize) {
	analyser_ = new VelocityField(info, sele1, voxelSize);
	break;
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'velocityField' requires the follwing specified\n"
		"\tSele1, VoxelSize\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerVelocityZ:
    default:
      break;
    }

    
    queryTime_ = analyzerParams->getQueryTime();
    std::cerr << "queryTime = " << queryTime_;
    

  }
    
  Analyzer::~Analyzer() {
    if (!doAnalyzer_) return;
#ifdef IS_MPI
    if (worldRank == 0) {
#endif

      
#ifdef IS_MPI
    }
#endif

  }

  
  void Analyzer::doAnalyzer() {
    if (!doAnalyzer_) return;
    // Here we call the StaticProps module on the single frame.
    // We may need to pass information here so we can pass it to
    //   writeOutputFile?
    // processFrame shouldn't need an int, check this.
    int itemp = 1;
    analyser_->processFrame(itemp);
  }
  
  
 
  void Analyzer::getStarted() {
    if (!doAnalyzer_) return;
    doAnalyzer();
    writeOutput();
  }
  
  void Analyzer::parseOutputFileFormat(const std::string& format) {
    if (!doAnalyzer_) return;
  }

  void Analyzer::writeOutput() {
    if (!doAnalyzer_) return;
    // Here we should call the writeOutput for the StaticProps modules.
    analyser_->writeOutput();
  }
   
    
}

