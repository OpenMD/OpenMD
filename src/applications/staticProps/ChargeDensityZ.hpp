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

#ifndef APPLICATIONS_STATICPROPS_CHARGEDENSITYZ_HPP
#define APPLICATIONS_STATICPROPS_CHARGEDENSITYZ_HPP
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "utils/ElementsTable.hpp"
#include "brains/Thermo.hpp"
#include "applications/staticProps/StaticAnalyser.hpp"
#include<set>
#include<map>

namespace OpenMD {

  class ChargeDensityZ : public StaticAnalyser{
  public:
    ChargeDensityZ(SimInfo* info, const std::string& filename,
		   const std::string& sele, int nzbins, RealType vRadius,std::string atomName = "Au", bool xyzGen=false, int axis=2);
    virtual void process();



  private:

    virtual void writeDensity();
    virtual void generateXYZForLastFrame();

    Snapshot* currentSnapshot_;
    int nProcessed_;
    std::string selectionScript_;
    SelectionEvaluator evaluator_;
    SelectionManager seleMan_;
    Thermo thermo_;

    std::vector<RealType> zBox_;
    std::vector<RealType> densityZAverageAllFrame_;
    std::vector<RealType> averageDensityZ_;
    std::vector<RealType> flucDensityZAverageAllFrame_;
    std::vector<RealType> densityFlucZAverageAllFrame_;
    std::vector<RealType> absDensityFlucZAverageAllFrame_;
    std::vector<RealType> densityFlucZAverageFirstFrame_;
    std::vector<RealType> absDensityFlucZAverageFirstFrame_;

    int axis_,x_,y_;
    RealType vRadius_;
    std::string fileName_;
    std::string atomFlucCharge_;
    bool genXYZ_;

    Mat3x3d hmat_;

    std::map<std::string,RealType> vander_waals_r;
    std::map<std::string, RealType> averageChargeForEachType_;
    std::map<std::string, int> SDCount_;
    std::string axisLabel_;

    std::map<int,RealType> averageChargeUsingGlobalIndex_;
    std::map<int,std::vector<RealType> > totalChargeUsingGlobalIndex_;
    std::map<int,RealType > totalChargeFluctationsUsingGlobalIndex_;
    std::map<int,std::vector<RealType> > zPosUsingGlobalIndex_;
    std::map<int,std::vector<RealType> > xPosUsingGlobalIndex_;
    std::map<int,std::vector<RealType> > yPosUsingGlobalIndex_;
    std::map<int,int> countUsingGlobalIndex_;
    std::map<int,RealType> vanderRUsingGlobalIndex_;
    std::map<int, std::string> atomNameGlobalIndex_;
  };

}

#endif
