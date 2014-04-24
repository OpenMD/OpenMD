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
 *  @version $Id$
 *
 */

#ifndef APPLICATIONS_STATICPROPS_BOPOFR_HPP
#define APPLICATIONS_STATICPROPS_BOPOFR_HPP
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "applications/staticProps/StaticAnalyser.hpp"
#include "math/Vector3.hpp"
#include "math/SphericalHarmonic.hpp"

namespace OpenMD {

  class BOPofR : public StaticAnalyser{
  public:
    BOPofR(SimInfo* info, const std::string& filename, 
           const std::string& sele, double rCut, int nbins, RealType len);
    
    virtual ~BOPofR();
    virtual void process();
    
  protected:
    virtual void initializeHistogram();
    virtual void collectHistogram(std::vector<RealType> q, 
                                  std::vector<ComplexType> what, 
                                  RealType distCOM) = 0;
    virtual void writeOrderParameter() = 0;

    Snapshot* currentSnapshot_;
    std::string selectionScript_;
    SelectionManager seleMan_;    
    SelectionEvaluator evaluator_;           
            
    RealType rCut_;
    static const int lMax_ = 6;
    int frameCounter_;
    int nBins_;
    RealType len_;
    RealType deltaR_;
    
    std::map<std::pair<int,int>,int> m2Min;
    std::map<std::pair<int,int>,int> m2Max;
    std::map<std::pair<int,int>,std::vector<RealType> > w3j;
   
    std::vector<int> RCount_;
    std::vector<int> WofR_;
    std::vector<int> QofR_;
  };


  class IcosahedralOfR : public BOPofR {
  public:
    IcosahedralOfR(SimInfo* info, const std::string& filename, 
                   const std::string& sele, double rCut, int nbins, 
                   RealType len);

    virtual void collectHistogram(std::vector<RealType> q, 
                                  std::vector<ComplexType> what, 
                                  RealType distCOM);
    virtual void writeOrderParameter();
  };

  class FCCOfR : public BOPofR {
  public:
    FCCOfR(SimInfo* info, const std::string& filename, 
           const std::string& sele, double rCut, int nbins, RealType len);

    virtual void collectHistogram(std::vector<RealType> q, 
                                  std::vector<ComplexType> what, 
                                  RealType distCOM);
    virtual void writeOrderParameter();
  };  
}

#endif

