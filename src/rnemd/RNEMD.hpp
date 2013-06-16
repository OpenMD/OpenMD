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
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
/**
 * @file RNEMD.hpp
 * @author gezelter
 * @date 03/13/2009
 * @time 15:56pm
 * @version 1.0
 */

#ifndef INTEGRATORS_RNEMD_HPP
#define INTEGRATORS_RNEMD_HPP
#include "brains/SimInfo.hpp"
#include "math/RandNumGen.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include <iostream>

using namespace std;
namespace OpenMD {

  class RNEMD {
  public:
    RNEMD(SimInfo* info);
    ~RNEMD();
    
    void doRNEMD();
    void doSwap(SelectionManager& smanA, SelectionManager& smanB);
    void doNIVS(SelectionManager& smanA, SelectionManager& smanB);
    void doVSS(SelectionManager& smanA, SelectionManager& smanB);
    RealType getDividingArea();
    void collectData();
    void getStarted();
    void parseOutputFileFormat(const std::string& format);
    void writeOutputFile();
    void writeReal(int index, unsigned int bin);
    void writeVector(int index, unsigned int bin);
    void writeRealStdDev(int index, unsigned int bin);
    void writeVectorStdDev(int index, unsigned int bin);

  private:

    enum RNEMDMethod {
      rnemdSwap,
      rnemdNIVS,
      rnemdVSS,
      rnemdUnkownMethod
    };
    enum RNEMDFluxType {
      rnemdKE,       // translational kinetic energy flux
      rnemdRotKE,    // rotational kinetic energy flux
      rnemdFullKE,   // full kinetic energy flux
      rnemdPx,       // flux of momentum along x axis 
      rnemdPy,       // flux of momentum along y axis 
      rnemdPz,       // flux of momentum along z axis 
      rnemdPvector,  // flux of momentum vector
      rnemdLx,       // flux of angular momentum along x axis 
      rnemdLy,       // flux of angular momentum along y axis 
      rnemdLz,       // flux of angular momentum along z axis 
      rnemdLvector,  // flux of angular momentum vector
      rnemdKePx,     // flux of translational KE and x-momentum
      rnemdKePy,     // flux of translational KE and y-momentum
      rnemdKePvector, // full combo flying platter
      rnemdKeLx,     // flux of translational KE and x-angular momentum
      rnemdKeLy,     // flux of translational KE and y-angular momentum
      rnemdKeLz,     // flux of translational KE and z-angular momentum
      rnemdKeLvector, // full combo spinning platter
      rnemdUnknownFluxType
    };

    enum OutputFields {
      BEGININDEX = 0,
      Z = BEGININDEX,
      R,
      TEMPERATURE,
      VELOCITY,
      ANGULARVELOCITY,
      DENSITY,
      ENDINDEX 
    };

    struct OutputData {
      string title;
      string units;
      string dataType;
      vector<BaseAccumulator*> accumulator;
    };

    typedef bitset<ENDINDEX-BEGININDEX> OutputBitSet;
    typedef map<string, OutputFields> OutputMapType;
    
    SimInfo* info_;

    map<string, RNEMDMethod> stringToMethod_;
    map<string, RNEMDFluxType> stringToFluxType_;
    RNEMDMethod rnemdMethod_;
    RNEMDFluxType rnemdFluxType_;

    // object selection for specifying a particular species:
    string rnemdObjectSelection_;
    SelectionEvaluator evaluator_;
    SelectionManager seleMan_;

    // Geometric selections for the two regions for the exchange:
    string selectionA_;
    SelectionEvaluator evaluatorA_;
    SelectionManager seleManA_;
    string selectionB_;
    SelectionEvaluator evaluatorB_;
    SelectionManager seleManB_;
    SelectionManager commonA_;
    SelectionManager commonB_;
    bool hasSelectionA_;
    bool hasSelectionB_;                      
    bool hasSphereBRadius_;

    bool usePeriodicBoundaryConditions_;
    bool hasDividingArea_;
    RealType dividingArea_;

    int nBins_;
    RealType binWidth_;
    RealType slabWidth_;
    RealType slabACenter_;
    RealType slabBCenter_;
    RealType sphereARadius_;
    RealType sphereBRadius_;

    Vector3d coordinateOrigin_;

    RealType kineticFlux_;        // target or desired *flux*
    Vector3d momentumFluxVector_; // target or desired *flux*
    Vector3d angularMomentumFluxVector_; // target or desired *flux*

    RealType kineticTarget_;     // target or desired one-time exchange energy
    Vector3d momentumTarget_;    // target or desired one-time exchange momentum
    Vector3d angularMomentumTarget_; // target or desired one-time
                                     // exchange angular momentum

    RealType kineticExchange_;    // actual exchange energy (running total)
    Vector3d momentumExchange_;   // actual exchange momentum (running total)
    Vector3d angularMomentumExchange_; // actual exchange momentum
                                       // (running total)
    RealType exchangeTime_;

    RealType targetJzpz2_;

    unsigned int trialCount_;
    unsigned int failTrialCount_;
    unsigned int failRootCount_;

    string rnemdFileName_;
    ofstream rnemdFile_;

    RealType runTime_, statusTime_;

    vector<OutputData> data_;
    OutputBitSet outputMask_;
    OutputMapType outputMap_;
    Accumulator* areaAccumulator_;
    bool doRNEMD_;
    bool hasData_;

  };
}
#endif //INTEGRATORS_RNEMD_HPP
