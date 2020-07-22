/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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

/*
* Computes the Cos\theta distribution along preferred axis for the selected atom
*/

#include <algorithm>
#include <fstream>
#include "applications/staticProps/OrderParameterProbZ.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "brains/Thermo.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"


namespace OpenMD {

  OrderParameterProbZ::OrderParameterProbZ(SimInfo* info, const std::string& filename,const std::string& sele,
    const RealType dipoleX, const RealType dipoleY, const RealType dipoleZ,
    int nbins, int axis)
    :StaticAnalyser(info, filename, nbins),
    selectionScript_(sele),evaluator_(info), seleMan_(info), thermo_(info), nbins_(nbins), axis_(axis) {

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // fixed number of bins
    Count_.resize(nbins);
    std::fill(Count_.begin(), Count_.end(), 0);

    switch(axis_) {
    case 0:
      axisLabel_ = "x";
      refAxis_ = Vector3d(1,0,0);
      break;
    case 1:
      axisLabel_ = "y";
      refAxis_ = Vector3d(0,1,0);
      break;
    case 2:
    default:
      axisLabel_ = "z";
      refAxis_ = Vector3d(0,0,1);
      break;
    }


    dipoleVector_ = Vector3d(dipoleX, dipoleY, dipoleZ);
    dipoleVector_.normalize();

    setOutputName(getPrefix(filename) + ".OrderProb");
  }

  void OrderParameterProbZ::process() {
    StuntDouble* sd;
    int ii;
    RealType orderMin = -1.0;
    RealType orderMax = 1.0;
    RealType deltaOrder = (orderMax - orderMin)/nbins_;

    bool usePeriodicBoundaryConditions_ =
      info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    totalCount_ = 0;
    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();


      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      //wrap the stuntdoubles into a cell and find order parameter

      for (sd = seleMan_.beginSelected(ii); sd != NULL; sd = seleMan_.nextSelected(ii)) {
      Vector3d pos = sd->getPos();
      if (usePeriodicBoundaryConditions_)
          currentSnapshot_->wrapVector(pos);
      sd->setPos(pos);
    }
    SquareMatrix3<RealType> rotMat;
    Vector3d rotatedDipoleVector;
    RealType ctheta;
    for (sd = seleMan_.beginSelected(ii); sd != NULL; sd = seleMan_.nextSelected(ii)) {
      if (sd->isDirectional() || sd->isRigidBody()) {
        rotMat = sd->getA();
        rotatedDipoleVector = rotMat * dipoleVector_;
        rotatedDipoleVector.normalize();
        ctheta = dot(rotatedDipoleVector, refAxis_);
        int index = int( (ctheta - orderMin)/deltaOrder );
        Count_[index]++;
        totalCount_++;
      }

  }

  }


    writeOrderCount();

}

  void OrderParameterProbZ::writeOrderCount() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#Order count probablity "<<"\n";
      rdfStream << "#selection: (" << selectionScript_ << ")\n";
      rdfStream << "# Prefered Axis:" << axisLabel_ << "\n##Order\tProbOrderCount\n";
      for (unsigned int i = 0; i < Count_.size(); ++i) {
        RealType order = i * (2.0/Count_.size());
        RealType prop;
        if (totalCount_ == 0) prop = Count_[i];
        else prop = Count_[i]/totalCount_;
        rdfStream << order << "\t"
                  << prop <<"\n";

      }

    } else {

      sprintf(painCave.errMsg, "OrderProb: unable to open %s\n",
	      outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();

  }
}
