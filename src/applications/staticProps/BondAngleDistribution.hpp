/*
 * Copyright (c) 2007 The University of Notre Dame. All Rights Reserved.
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
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [4] , Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). *
 *  Created by J. Daniel Gezelter on 07/27/07
 *  @author  J. Daniel Gezelter 
 *  @version $Id$
 *
 */

#ifndef APPLICATIONS_STATICPROPS_BONDANGLEDISTRIBUTION_HPP
#define APPLICATIONS_STATICPROPS_BONDANGLEDISTRIBUTION_HPP
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "applications/staticProps/StaticAnalyser.hpp"


namespace OpenMD {

  /**
   * @class BondAngleDistribution
   * @brief BondAngleDistribution
   *
   *  Comptes bond angle distribution for nearest neighbors.
   *BondAngleDistribution
   */
  class BondAngleDistribution : public StaticAnalyser{
  public:
    BondAngleDistribution(SimInfo* info, const std::string& filename, 
                       const std::string& sele, double rCut, int nbins);
    
    virtual ~BondAngleDistribution();
    virtual void process();
    
  private:
    virtual void initializeHistogram();
       
    void writeBondAngleDistribution();

    Snapshot* currentSnapshot_;
    std::string selectionScript_;
    SelectionManager seleMan_;    
    SelectionEvaluator evaluator_;           
            
    RealType rCut_;
    int frameCounter_;
    int nBins_;
    

    RealType deltaTheta_;
    std::vector<int> ThetaCount_;
    std::vector<int> histogram_;
    int nTotBonds_;
   

  };
}

#endif

