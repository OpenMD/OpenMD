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

#include "analyzer/Analyzer.hpp"
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

#include "applications/staticProps/StaticAnalyser.hpp"
#include "applications/staticProps/GofR.hpp"
#include "applications/staticProps/GofZ.hpp"
#include "applications/staticProps/GofRZ.hpp"
#include "applications/staticProps/GofRAngle.hpp"
#include "applications/staticProps/GofAngle2.hpp"
#include "applications/staticProps/GofRAngle2.hpp"
#include "applications/staticProps/GofXyz.hpp"
#include "applications/staticProps/TwoDGofR.hpp"
#include "applications/staticProps/P2OrderParameter.hpp"
#include "applications/staticProps/BondOrderParameter.hpp"
#include "applications/staticProps/BOPofR.hpp"
#include "applications/staticProps/RippleOP.hpp"
#include "applications/staticProps/SCDOrderParameter.hpp"
#include "applications/staticProps/DensityPlot.hpp"
#include "applications/staticProps/ObjectCount.hpp"
#include "applications/staticProps/RhoZ.hpp"
#include "applications/staticProps/PipeDensity.hpp"
#include "applications/staticProps/pAngle.hpp"
#include "applications/staticProps/BondAngleDistribution.hpp"
#include "applications/staticProps/NanoVolume.hpp"
#include "applications/staticProps/NanoLength.hpp"
#if defined(HAVE_FFTW_H) || defined(HAVE_DFFTW_H) || defined(HAVE_FFTW3_H)
#include "applications/staticProps/Hxy.hpp"
#endif
#include "applications/staticProps/RhoR.hpp"
#include "applications/staticProps/AngleR.hpp"
#include "applications/staticProps/TetrahedralityParam.hpp"
#include "applications/staticProps/TetrahedralityParamZ.hpp"
#include "applications/staticProps/TetrahedralityParamXYZ.hpp"
#include "applications/staticProps/TetrahedralityParamDens.hpp"
#include "applications/staticProps/RNEMDStats.hpp"
#include "applications/staticProps/NitrileFrequencyMap.hpp"
#include "applications/staticProps/MultipoleSum.hpp"
#include "applications/staticProps/SurfaceDiffusion.hpp"
#include "applications/staticProps/CoordinationNumber.hpp"
#include "applications/staticProps/HBondGeometric.hpp"
#include "applications/staticProps/PotDiff.hpp"
#include "applications/staticProps/TetrahedralityHBMatrix.hpp"
#include "applications/staticProps/Kirkwood.hpp"
#include "applications/staticProps/Field.hpp"
#include "applications/staticProps/VelocityZ.hpp"


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
      RealType step = analyzerParams->getStep();
    if (hasNBins)
      RealType nbins = analyzerParams->getNBins();
    if (hasNBinsX)
      RealType nbinsX = analyzerParams->getNBinsX();
    if (hasNBinsY)
      RealType nbinsY = analyzerParams->getNBinsY();
    if (hasNBinsZ)
      RealType nbinsZ = analyzerParams->getNBinsZ();
    if (hasNRBins)
      RealType nrbins = analyzerParams->getNRBins();
    if (hasNAngleBins)
      RealType nAngleBins = analyzerParams->getNAngleBins();
    if (hasRCut)
      RealType rcut = analyzerParams->getRCut();
    if (hasOOCut)
      RealType ooCut = analyzerParams->getOOCut();
    if (hasThetaCut)
      RealType thetaCut = analyzerParams->getThetaCut();
    if (hasOHCut)
      RealType OHCut = analyzerParams->getOHCut();
    if (hasDz)
      RealType dz = analyzerParams->getDz();
    if (hasLength)
      RealType length = analyzerParams->getLength();
    if (hasZLength)
      RealType zLength = analyzerParams->getZLength();
    if (hasZOffSet)
      RealType zOffSet = analyzerParams->getZOffSet();
    if (hasSele1)
      std::string sele1 = analyzerParams->getSele1();
    if (hasSele2)
      std::string sele2 = analyzerParams->getSele2();
    if (hasSele3)
      std::string sele3 = analyzerParams->getSele3();
    if (hasComSele)
      std::string comSele = analyzerParams->getComSele();
    if (hasSeleOffSet)
      RealType seleOffSet = analyzerParams->getSeleOffSet();
    if (hasSeleOffSet2)
      RealType seleOffSet2 = analyzerParams->getSeleOffSet2();
    if (hasMolName)
      std::string molName = analyzerParams->getMolName();
    if (hasBegin)
      RealType begin = analyzerParams->getBegin();
    if (hasEnd)
      RealType end = analyzerParams->getEnd();
    if (hasRadius)
      RealType radius = analyzerParams->getRadius();
    if (hasVoxelSize)
      RealType voxelSize = analyzerParams->getVoxelSize();
    if (hasGaussWidth)
      RealType gaussWidth = analyzerParams->getGaussWidth();
    if (hasPrivilegedAxis)
      std::string privilegedAxis = analyzerParams->getPrivilegedAxis();
    if (hasPrivilegedAxis2)
      std::string privilegedAxis2 = analyzerParams->getPrivilegedAxis2();

    
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

    
    StaticAnalyser* analyser;
    
    switch(analyzerMethod_) {
    case analyzerBo:
      if (hasSele1 && hasSele2 && hasLength && hasNRBins) {
	//analyser= new GofR(info, sele1, sele2, maxLen, nrbins);
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'gofr' requires the following specified\n"
		"\tSele1, Sele2, Length, and NRBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerIor:
      if (hasSele1 && hasSele2 && hasLength && hasNBins && hasPrivilegedAxis) {
	//analyser = new GofZ(info, sele1, sele2, length, nbins, privilegedAxis);
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'gofz' requires the following specified\n"
		"\tSele1, Sele2, Length, NBins, privilegedAxis\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerFor:
      if (hasSele1 && hasRCut && hasNRBins && hasLength && hasPrivilegedAxis) {
	//analyser = new FCCofR(info, sele1, rcut, nrbins, length);
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'for' requires the follwing specified\n"
		"\tSele1, rCut, nRBins, Length\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerBad:
      if (hasSele1 && hasRCut && hasNBins) {
	// analyser = new BondAngleDistribution(info, sele1, rcut, nbins);
      } else {
	sprintf(painCave.errMsg,
		"Analyzer: 'bad' requires the follwing specified\n"
		"\tSele1, rCut, nBins\n");
	painCave.isFatal = 1;
	painCave.severity = OPENMD_ERROR;
	simError();
      }
    case analyzerCount:
      if (hasSele1) {
	//analyser = new ObjectCount(info, sele1);
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
	// analyser = new GofR(info, sele1, sele2, length, nrbins);
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
	// analyser = new GofZ(info, sele1, sele2, length, nbins, privilegedAxis);
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
	// analyser = new GofRTheta(info, sele1, sele2, sele3, length, nrbins, nanglebins);
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
	// analyser = new GofROmega(info, sele1, sele2, sele3, length, nrbins, nanglebins);
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
	// analyser = new GofRZ(info, sele1, sele2, length, maxLength, nrbins, nbinsZ, privilegedAxis);
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
	// anlyser = new GofAngle2(info, sele1, sele2, sele3, nagnlebins);
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
	// analyser = new GofRAngle2(info, sele1, sele2, sele3, maxLength, nrbins, nanglebins);
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
	// analyser = new GofXyz(info, sele1, sele2, maxLength, nbins);
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
	// analyser = new TwoDGofR(indo, sele1, sele2, maxLength, dz, nrbins);
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
	  // analyser = new P2OrderParameter(info, sele1, sele2);
	}
	else {
	  if (hasSeleOffSet) {
	    // analyser = new P2OrderParameter(info, sele1, seleOffSet);
	  }
	  else {
	    // analyser = new P2OrderParameter(info, sele1);
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
	// analyser = new SCDOrderParameter(info, sele1, sele2, sele3);
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
	// analyser = new DensityPlot(info, sele1, sele2, maxLength, nbins);
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
	// analyser = new RhoZ(info, sele1, nbins, privilegedAxis);
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
	  // analyser = new pAngle(info, sele1, sele2, nbins);
	} else {
	  if (hasSeleOffSet) {
	    if (hasSeleOffSet2) {
	      // analyser = new pAngle(info, sele1, seleOffSet, seleOffSet2, nbins);
	    } else {
	      // analyser = new pAngle(info, sele1, seleOffSet, nbins);
	    }
	  }
	  // analyser = new pAngle(info, sele1, nbins);
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
	// analyser = new Hxy(info, sele1, nbinsX, nbinsY, nbinsZ, nbins);
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
	// analyser = new RhoR(info, sele1, maxLength, nrbins, radius);
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
	// analyser = new AngleR(info, sele1, maxLength, nrbins);
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
	// analyser = new NanoVolume(info, sele1);
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
	// analyser = new NanoLength(info, sele1);
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
	// analyser = new TetrahedralityParam(info, sele1, rcut, nbins);
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
	// analyser = new TetrahedralityParamZ(info, sele1, sele2, rcut, nbins, privilegedAxis);
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
	// analyser = new TetrahedralityParamDens(info, sele1, sele2, rcut, nbins);
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
	// analyser = new TetrahedralityParamXYZ(info, sele1, sele2, rcut, voxelSize, gaussWidth);
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
	// analyser = new RNEMDZ(info, sele1, nbins, privilegedAxis);
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
	// analyser = new RNEMDR(info, sele1, nrbins);
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
	// analyser = new RNEMDRTheta(info, sele1, nrbins, nanglebins);
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
	// analyser = new NitrileFrequencyMap(info, sele1, nbins);
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
	// analyser = new MultipoleSum(info, sele1, maxLength, nbins);
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
	// analyser = new SurfaceDiffusion(info, sele1, maxLength);
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
	// analyser = new CoordinationNumber(info, sele1, sele2, rcut, nbins);
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
	// analyser = new SCN(info, sele1, sele2, rcut, nbins);
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
	// analyser = new GCN(info, sele1, sele2, rcut, nbins);
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
	// analyser = new HBondGeometric(info, sele1, sele2, rcut, thetaCut, nbins);
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
	// analyser = new PotDiff(info, sele1);
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
	// analyser = new TetrahedralityHBMatrix(info, sele1, rcut, ooCut, thetaCut, OHCut, nbins);
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
	// analyser = new Kirkwood(info, sele1, sele2, maxLength, nrbins);
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
	// analyser = new KirkwoodQ(info, sele1, sele2, maxLength, nrbins);
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
	// analyser = new DensityField(info, sele1, voxelSize);
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
	// analyser = new VelocityField(info, sele1, voxelSize);
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
  }
  
  
  void Analyzer::collectData() {
    if (!doAnalyzer_) return;
  }
  
 
  void Analyzer::getStarted() {
    if (!doAnalyzer_) return;
    collectData();
    writeOutputFile();
  }
  
  void Analyzer::parseOutputFileFormat(const std::string& format) {
    if (!doAnalyzer_) return;
  }

  void Analyzer::writeOutputFile() {
    if (!doAnalyzer_) return;
  }
   
    
}

