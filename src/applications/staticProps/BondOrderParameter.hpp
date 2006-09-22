/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
#ifndef APPLICATIONS_STATICPROPS_BONDORDERPARAMETER_HPP
#define APPLICATIONS_STATICPROPS_BONDORDERPARAMETER_HPP
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "applications/staticProps/StaticAnalyser.hpp"
#include "math/Vector3.hpp"
#include "math/Wigner3jm_interface.h"

namespace oopse {
  
  class BondOrderParameter : public StaticAnalyser{
  public:
    BondOrderParameter(SimInfo* info, const std::string& filename, 
                       const std::string& sele, double rCut, int lNumber, int nbins);

    virtual ~BondOrderParameter();
    virtual void process();

  private:
    
    void writeOrderParameter(RealType Q_l, RealType W_l_hat);
    virtual void initalizeHistogram();
    virtual void collectHistogram(RealType q_l, RealType w_l);

    Snapshot* currentSnapshot_;

    std::string selectionScript_;
    SelectionManager seleMan_;    
    SelectionEvaluator evaluator_;           
            
    RealType rCut_;
    int lNumber_;
    int mSize_;    
    int frameCounter_;

    RealType MinQ_;
    RealType MaxQ_;
    RealType deltaQ_;
    RealType sumQ_;
    RealType sumQ2_;
    int Qcount_;
    std::vector<int> Q_histogram_;

    RealType MinW_;
    RealType MaxW_;
    RealType deltaW_;
    RealType sumW_;
    RealType sumW2_;
    int Wcount_;
    std::vector<int> W_histogram_;
  };
}

#endif

