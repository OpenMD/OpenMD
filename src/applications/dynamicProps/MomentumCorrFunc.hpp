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

 /* Uses the Helfand-moment method for calculating thermal
  * conductivity using the relation kappa = (N,V)lim(t)->inf 1/(2*k_B*T^2*V*t) <[G_K(t)-G_K(0)]^2>
  * where G_K is the Helfand moment for thermal conductivity definded as
  * G_K(t) = sum_{a=1}{^N} x_a(E_a-<E_a>) and E_a is defined to be
  *	E_a = p_2^2/(2*m)+1/2 sum_{b.ne.a} u(r_ab) where p is momentum and u is pot energy for the 
  * particle pair a-b. This routine calculates E_a, <E_a> and does the correlation
  * <[G_K(t)-G_K(0)]^2>.
  * See Viscardy et al. JCP 126, 184513 (2007)
  */


#ifndef APPLICATIONS_DYNAMICPROPS_MOMENTUMCORRFUNC_HPP
#define APPLICATIONS_DYNAMICPROPS_MOMENTUMCORRFUNC_HPP

#include "applications/dynamicProps/FrameTimeCorrFunc.hpp"

namespace OpenMD {

  class MomentumCorrFunc : public FrameTimeCorrFunc {
  public:
    MomentumCorrFunc(SimInfo* info, const std::string& filename, const std::string& sele1, const std::string& sele2, long long int memSize);   
        
  private:
    virtual void correlateFrames(int frame1, int frame2);
    virtual RealType calcCorrVal(int frame1, int frame2) { return 0.0; }
    virtual void writeCorrelate();

  protected:
    virtual void preCorrelate();
    virtual void postCorrelate();
    std::vector<Mat3x3d> histogram_;
    
  };

}
#endif //MOMENTUMCORRFUNC


