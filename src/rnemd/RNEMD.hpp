/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#ifndef OPENMD_RNEMD_RNEMD_HPP
#define OPENMD_RNEMD_RNEMD_HPP

#include <bitset>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "brains/SimInfo.hpp"
#include "brains/Snapshot.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "types/AtomType.hpp"
#include "utils/Accumulator.hpp"
#include "utils/StaticAccumulator.hpp"

namespace OpenMD {
  namespace RNEMD {

    class RNEMD {
    public:
      RNEMD(SimInfo* info);
      virtual ~RNEMD();

      void getStarted();
      void doRNEMD();
      void collectData();
      void writeOutputFile();

    protected:
      enum RNEMDFluxType {
        rnemdKE,         // translational kinetic energy flux
        rnemdRotKE,      // rotational kinetic energy flux
        rnemdFullKE,     // full kinetic energy flux
        rnemdPx,         // flux of momentum along x axis
        rnemdPy,         // flux of momentum along y axis
        rnemdPz,         // flux of momentum along z axis
        rnemdPvector,    // flux of momentum vector
        rnemdLx,         // flux of angular momentum along x axis
        rnemdLy,         // flux of angular momentum along y axis
        rnemdLz,         // flux of angular momentum along z axis
        rnemdLvector,    // flux of angular momentum vector
        rnemdKePx,       // flux of translational KE and x-momentum
        rnemdKePy,       // flux of translational KE and y-momentum
        rnemdKePvector,  // full combo flying platter
        rnemdKeLx,       // flux of translational KE and x-angular momentum
        rnemdKeLy,       // flux of translational KE and y-angular momentum
        rnemdKeLz,       // flux of translational KE and z-angular momentum
        rnemdKeLvector,  // full combo spinning platter
        rnemdUnknownFluxType
      };

      Snapshot* currentSnap_ {nullptr};

      RNEMDFluxType rnemdFluxType_;

      std::string rnemdFluxTypeLabel_;
      std::string rnemdMethodLabel_;

      bool doRNEMD_ {false};
      bool usePeriodicBoundaryConditions_ {false};

      Vector3d coordinateOrigin_;

      RealType kineticTarget_ {
          0.0};  // target or desired one-time exchange energy
      Vector3d momentumTarget_ {
          V3Zero};  // target or desired one-time exchange momentum
      Vector3d angularMomentumTarget_ {
          V3Zero};  // target or desired one-time exchange angular momentum

      RealType kineticExchange_ {
          0.0};  // actual exchange energy (running total)
      Vector3d momentumExchange_ {
          V3Zero};  // actual exchange momentum (running total)
      Vector3d angularMomentumExchange_ {
          V3Zero};  // actual exchange momentum (running total)

      unsigned int trialCount_ {0};
      unsigned int failTrialCount_ {0};
      unsigned int failRootCount_ {0};

      void setKineticFlux(RealType kineticFlux);
      void setMomentumFluxVector(
          const std::vector<RealType>& momentumFluxVector);
      void setAngularMomentumFluxVector(
          const std::vector<RealType>& angularMomentumFluxVector);

    private:
      enum OutputFields {
        BEGININDEX = 0,
        Z          = BEGININDEX,
        R,
        TEMPERATURE,
        VELOCITY,
        ANGULARVELOCITY,
        DENSITY,
        ACTIVITY,
        ELECTRICFIELD,
        ELECTROSTATICPOTENTIAL,
        ENDINDEX
      };

      enum RNEMDPrivilegedAxis { rnemdX = 0, rnemdY = 1, rnemdZ = 2 };

      struct OutputData {
        std::string title;
        std::string units;
        std::string dataType;
        std::vector<BaseAccumulator*> accumulator;
        std::vector<std::vector<BaseAccumulator*>> accumulatorArray2d;
      };

      using OutputBitSet  = std::bitset<ENDINDEX - BEGININDEX>;
      using OutputMapType = std::map<std::string, OutputFields>;

      SimInfo* info_ {nullptr};
      Mat3x3d hmat_;

      RNEMDPrivilegedAxis rnemdPrivilegedAxis_;
      std::string rnemdAxisLabel_;

      // Object selection for specifying a particular species:
      std::string rnemdObjectSelection_;
      SelectionEvaluator evaluator_;
      SelectionManager seleMan_;

      // Geometric selections for the two regions for the exchange:
      std::string selectionA_;
      std::string selectionB_;
      SelectionEvaluator evaluatorA_;
      SelectionEvaluator evaluatorB_;
      SelectionManager seleManA_;
      SelectionManager seleManB_;
      SelectionManager commonA_;
      SelectionManager commonB_;
      bool hasSelectionA_ {false};
      bool hasSelectionB_ {false};

      // Output selection for collecting data about a particular species:
      std::string outputSelection_;
      SelectionEvaluator outputEvaluator_;
      SelectionManager outputSeleMan_;

      unsigned int nBins_ {0};

      RealType binWidth_;
      RealType slabWidth_;

      RealType slabACenter_;
      RealType slabBCenter_;
      RealType sphereARadius_;
      RealType sphereBRadius_;
      RealType areaA_;
      RealType areaB_;
      RealType volumeA_;
      RealType volumeB_;

      bool hasData_ {false};
      bool hasDividingArea_ {false};

      RealType dividingArea_;

      RealType kineticFlux_ {0.0};                   // target or desired *flux*
      Vector3d momentumFluxVector_ {V3Zero};         // target or desired *flux*
      Vector3d angularMomentumFluxVector_ {V3Zero};  // target or desired *flux*

      RealType exchangeTime_ {}, runTime_ {}, statusTime_ {};

      std::string rnemdFileName_;
      std::ofstream rnemdFile_ {};

      OutputBitSet outputMask_;
      OutputMapType outputMap_;
      std::vector<OutputData> data_;
      std::vector<AtomType*> outputTypes_;

      Utils::RealAccumulator areaAccumulator_ {};

      void parseOutputFileFormat(const std::string& format);
      std::string setSelection(RealType& slabCenter);
      RealType getDefaultDividingArea();
      int getBin(Vector3d pos);
      void writeReal(int index, unsigned int bin);
      void writeVector(int index, unsigned int bin);
      void writeArray(int index, unsigned int bin);
      void writeRealErrorBars(int index, unsigned int bin);
      void writeVectorErrorBars(int index, unsigned int bin);
      void writeArrayErrorBars(int index, unsigned int bin);

      virtual void doRNEMDImpl(SelectionManager& smanA,
                               SelectionManager& smanB) = 0;
    };
  }  // namespace RNEMD
}  // namespace OpenMD

#endif  // OPENMD_RNEMD_RNEMD_HPP
