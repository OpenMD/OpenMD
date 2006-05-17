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
 *
 *
 *  hxy.hpp
 *  OOPSE-2.0
 *
 *  Created by Xiuquan Sun on 05/09/06.
 *  @author  Xiuquan Sun
 *  @version $Id: Hxy.hpp,v 1.4 2006-05-17 21:51:42 tim Exp $
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
namespace oopse {
  
  class Hxy : public StaticAnalyser {
    
  public:
    Hxy(SimInfo* info, const std::string& filename, const std::string& sele, int nbins_x, int nbins_y, int nrbins);
    
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
    
    RealType nBinsX_;
    RealType nBinsY_;
    int nbins_; 
    RealType dfreq;
    
    std::vector<RealType> gridZ_;
    std::vector<int> gridsample_;
    std::vector< std::vector<RealType> > bin;
    std::vector< std::vector<int> > samples;
    std::vector<RealType> sum_bin, sum_bin_sq, avg_bin, avg_bin_sq;
    std::vector<RealType> errbin_sum, errbin_sum_sq, errbin, errbin_sq;
  };
  
}
#endif



