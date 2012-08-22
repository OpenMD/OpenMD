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
 *  Created by Xiuquan Sun on 05/09/06.
 *  @author  Xiuquan Sun
 *  @version $Id$
 *
 */
#ifndef APPLICATIONS_STATICPROPS_HXY_HPP
#define APPLICATIONS_STATICPROPS_HXY_HPP

#include "config.h"

#if defined(HAVE_FFTW_H) || defined(HAVE_DFFTW_H) || defined(HAVE_FFTW3_H)

#ifdef HAVE_FFTW_H
#include <fftw.h>
#endif

#ifdef HAVE_DFFTW_H
#include <dfftw.h>
#endif

#ifdef HAVE_FFTW3_H
#include <fftw3.h>
#define c_re(c) ((c)[0])
#define c_im(c) ((c)[1])
#endif
#endif

#include "applications/staticProps/RadialDistrFunc.hpp"
namespace OpenMD {
  
  class Hxy : public StaticAnalyser {
    
  public:
    Hxy(SimInfo* info, const std::string& filename, const std::string& sele, int nbins_x, int nbins_y, int nrbins);
    virtual ~Hxy();
    int getnbins() {
      return nbins_; 
    }
    
    int getnBinsX() {
      return nBinsX_;
    }
    
    int getnBinsY() {
      return nBinsY_;
    }
    
    virtual void process();
    
  private:
    
    virtual void printSpectrum();
        
    Snapshot* currentSnapshot_;
    
    int nProcessed_;
    std::string selectionScript_;
    SelectionEvaluator evaluator_;
    SelectionManager seleMan_;
    
    int nBinsX_;
    int nBinsY_;
    int nbins_; 
    RealType dfreq;
    
    std::vector<RealType> gridZ_;
    std::vector<int> gridsample_;
    std::vector< std::vector<RealType> > bin;
    std::vector< std::vector<int> > samples;

    std::vector<RealType> sum_bin, sum_bin_sq, avg_bin, avg_bin_sq;
    std::vector<RealType> errbin_sum, errbin_sum_sq, errbin, errbin_sq;
    std::vector<RealType> mag, newmag;

  };
  
}
#endif



