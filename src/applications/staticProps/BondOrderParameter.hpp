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
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [4] , Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). *
 *  Created by J. Daniel Gezelter on 09/26/06
 *  @author  J. Daniel Gezelter 
 *  @version $Id$
 *
 */

#ifndef APPLICATIONS_STATICPROPS_BONDORDERPARAMETER_HPP
#define APPLICATIONS_STATICPROPS_BONDORDERPARAMETER_HPP
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "applications/staticProps/StaticAnalyser.hpp"
#include "math/Vector3.hpp"
#include "math/SphericalHarmonic.hpp"

namespace OpenMD {

  /**
   * @class BondOrderParameter
   * @brief Bond Order Parameter
   *
   * Computes orientational bond order parameters as outlined in:
   *
   *   "Bond-orientaional order in liquids and glasses," by
   *    P. J. Steinhart, D. R. Nelson, and M. Ronchetti, 
   *    Phys. Rev. B, 28, 784 (1983).
   *
   * A somewhat more useful reference which has formulae for these order
   * parameters for individual atoms is:
   *
   *   "Numerical calculation of the rate of crystal nucleation in a
   *    Lennard-Jones system at moderate undercooling," by
   *    Pieter Rein ten Wolde, Maria J. Ruiz-Montero, and Daan Frenkel,
   *    J. Chem. Phys. 104, pp. 9932-9947 (1996).
   *
   * Note that this version uses a single cutoff radius to decide 
   * membership in the list of neighbors, and does not have use a
   * distance-dependent weighting as used in the second reference above.
   *
   * The selection script can be utilized to look at specific types of
   * central atoms.  A dynamic selector can also be utilized.  By
   * default, this class computes the \f[ Q_{l} \f] and 
   * \f[ \hat{W}_{l} \f] parameters up to \f[ l = 12 \f].  The completed
   * configurational averages of these values as well as the
   * distributions of atomic \f[ q_{l} \f] and \f[ \hat{w}_{l} \f]
   * values are then placed in .boq and .bow files.
   */
  class BondOrderParameter : public StaticAnalyser{
  public:
    BondOrderParameter(SimInfo* info, const std::string& filename, 
                       const std::string& sele, double rCut, int nbins);
    
    virtual ~BondOrderParameter();
    virtual void process();
    
  private:
    virtual void initializeHistogram();
    virtual void collectHistogram(std::vector<RealType> q, 
                                  std::vector<ComplexType> what);    
    void writeOrderParameter(std::vector<RealType> Q, 
                             std::vector<ComplexType> What);

    Snapshot* currentSnapshot_;
    std::string selectionScript_;
    SelectionManager seleMan_;    
    SelectionEvaluator evaluator_;           
            
    RealType rCut_;
    static const int lMax_ = 12;
    int frameCounter_;
    int nBins_;
    
    std::map<std::pair<int,int>,int> m2Min;
    std::map<std::pair<int,int>,int> m2Max;
    std::map<std::pair<int,int>,std::vector<RealType> > w3j;
   
    RealType MinQ_;
    RealType MaxQ_;
    RealType deltaQ_;
    std::vector<int> Qcount_;
    std::map<std::pair<int,int>,int> Q_histogram_;

    RealType MinW_;
    RealType MaxW_;
    RealType deltaW_;
    std::vector<int> Wcount_;
    std::map<std::pair<int,int>,int> W_histogram_;
  };
}

#endif

