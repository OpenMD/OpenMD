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
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
#ifndef APPLICATIONS_STATICPROPS_RADIALDISTRFUNC_HPP
#define APPLICATIONS_STATICPROPS_RADIALDISTRFUNC_HPP

#include <string>
#include <vector>

#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "utils/NumericConstant.hpp"
#include "applications/staticProps/StaticAnalyser.hpp"

namespace OpenMD {

  /**
   * @class RadialDistrFunc
   * @brief Radial Distribution Function
   */
  class RadialDistrFunc : public StaticAnalyser {
  public:
    RadialDistrFunc(SimInfo* info, const std::string& filename, const std::string& sele1, const std::string& sele2);

    virtual ~RadialDistrFunc() {}
        
    void process();        


        
  protected:

    virtual void preProcess() {}
    virtual void postProcess() {}

    int getNPairs() { return nPairs_;}
        
    Snapshot* currentSnapshot_;

    std::string selectionScript1_;
    std::string selectionScript2_;
    int nProcessed_;
    SelectionManager seleMan1_;
    SelectionManager seleMan2_;
        
  private:

    virtual void initalizeHistogram() {}
    virtual void collectHistogram(StuntDouble* sd1, StuntDouble* sd2) =0;
    virtual void processHistogram() {}
    void processNonOverlapping(SelectionManager& sman1, SelectionManager& sman2);
    void processOverlapping(SelectionManager& sman);

    virtual void validateSelection1(SelectionManager& sman) {}
    virtual void validateSelection2(SelectionManager& sman) {}
    virtual void writeRdf() = 0;

        
    SelectionEvaluator evaluator1_;
    SelectionEvaluator evaluator2_;

    SelectionManager sele1_minus_common_;
    SelectionManager sele2_minus_common_;
    SelectionManager common_;        
    int nPairs_;
  };


}
#endif
