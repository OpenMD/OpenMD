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
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [4] , Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). *
 *  Created by J. Daniel Gezelter on 09/26/06
 *  @author  J. Daniel Gezelter 
 *  @version $Id: BondOrderParameter.hpp 1442 2010-05-10 17:28:26Z gezelter $
 *
 */

#ifndef ANALYSIS_HBONDGEOMETRIC_HPP
#define ANALYSIS_HBONDGEOMETRIC_HPP
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "analysis/NonSpatialStatistics.hpp"

namespace OpenMD {

  /**
   * @class HBondGeometric
   * @brief Hydrogen Bonding statistics using geometric criteria
   *
   * Computes hydrogen bonding statistics using structural information:
   *
   *   "Hydrogen-bond dynamics for the extended simple point-carge
   *   model of water" by F. W. Starr, J. K. Nielsen, and
   *   H. E. Stanley, Phys. Rev. E 62(1), 579-587 (2000).
   *
   */
  class HBondGeometric : public NonSpatialStatistics {
  public:

    HBondGeometric(SimInfo* info, 
                   const std::string& sele1, const std::string& sele2,
                   double rCut, double thetaCut,
                   int nbins=10);
    
    virtual ~HBondGeometric();
    virtual void processFrame(int frame);

  protected:
    OutputData* hBonds;
    OutputData* nAcceptor;
    OutputData* nDonor;
    OutputData* nSelected;
    
  private:
    virtual void processStuntDouble(StuntDouble* sd, int bin);  
    void writeOutput();

    std::string selectionScript1_;
    SelectionManager seleMan1_;    
    SelectionEvaluator evaluator1_;
    std::string selectionScript2_;
    SelectionManager seleMan2_;    
    SelectionEvaluator evaluator2_;
                
    RealType rCut_;
    RealType thetaCut_;
    int frameCounter_;
    int nBins_;
   
  };
}

#endif

