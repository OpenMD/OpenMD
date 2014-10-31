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
 *  Created by Joseph R. Michalka on 06/13/12.
 *  @author  Joseph R. Michalka
 *  @version $Id: SurfaceDisruption.hpp 1665 2011-11-22 20:38:56Z gezelter $
 *
 */
#ifndef APPLICATIONS_STATICPROPS_SURFACEDIFFUSION_HPP
#define APPLICATIONS_STATICPROPS_SURFACEDIFFUSION_HPP

#include <string>
#include <vector>
#include <fstream>
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "utils/NumericConstant.hpp"
#include "applications/staticProps/StaticAnalyser.hpp"

namespace OpenMD {
  
  class SurfaceDiffusion: public StaticAnalyser {
    
  public:
    SurfaceDiffusion(SimInfo* info, const std::string& filename, const std::string& sele, RealType len);
    virtual ~SurfaceDiffusion();
    
    virtual void process();
    
  private:
    
    double round(double r);
    void solventAccessible();
    void mobileAtomsFirst();
    void mobileAtomsLast();
    void mobileAtoms();
    void positionCorrelation();

    Snapshot* currentSnapshot_;
    
    int nProcessed_;
    std::string selectionScript_;
    SelectionEvaluator evaluator_;
    SelectionManager seleMan1_;
  
    string filename_;
    int bins_; 
    int selectionCount_;
    double singleMoveDistance_;
    int frames_;
    bool doSolvent_;


    //All positions of all frames of selected indices
    //positions_[0][i]
    // First particle at frame i
    std::vector< std::vector<Vector3d> > positions_;
    std::vector< std::vector<Vector3d> > positions2_;

    //mobility of particle i at time j
    // moBool[i][j]
    // moBool.resize(selectionCount_);
    // for(){
    //  moBool[i].resize(frames);
    // }
    std::vector< std::vector<bool> > moBool_;
    std::vector< std::vector<bool> > moBool2_;

    std::vector< std::vector<StuntDouble*> > gridSD_;
    std::vector<StuntDouble*> gridHighZ_;
    std::vector<StuntDouble*> gridLowZ_;
    std::vector<StuntDouble*> forIndex_;
    std::vector<Vector3d> firstPosition_;
    std::vector<Vector3d> lastPosition_;


    std::vector<int> SAIndices_;
    std::vector<int> mobileIndices_;
    std::vector<int> indices_;
    std::vector<int> count_;

    std::vector<RealType> xHist_;
    std::vector<RealType> yHist_;
    std::vector<RealType> rHist_;

    RealType minDistance_;
    RealType probe_;
  };
  
}
#endif



