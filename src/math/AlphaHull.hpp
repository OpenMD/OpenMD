/* Copyright (c) 2008, 2009 The University of Notre Dame. All Rights Reserved.
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
 *
 *  AlphaHull.hpp
 *
 *  Purpose: To calculate alphahull, hull volume using the QuickHull algorithm provided by QHull.
 *
 *  Created by Charles F. Vardeman II on 16 Dec 2009.
 *  @author  Charles F. Vardeman II
 *  @version $Id$
 *
 */

#ifndef MATH_ALPHAHULL_HPP_
#define MATH_ALPHAHULL_HPP_

#include "math/Vector3.hpp"
#include "config.h"
#include "math/Hull.hpp"
#include "math/Triangle.hpp"

#include <cassert>
#include <vector>
#include <string>


namespace OpenMD {
  class AlphaHull : public Hull {
  public:

    AlphaHull(double alpha);    
    virtual ~AlphaHull(){};

    void computeHull( std::vector<StuntDouble*> bodydoubles );

    /* Total area of Hull*/
    RealType getArea(){return area_;}

    /* Total Volume enclosed by Hull */
    RealType getVolume(){ return volume_; } 

    std::vector<Triangle> getMesh(){return Triangles_;}
    void printHull(const std::string& geomFileName);

  protected:
    RealType volume_;
    RealType area_;
    int dim_;
    double alpha_;
    const std::string options_;
    
  private:
    std::vector<Triangle> Triangles_;
  };
}
#endif /*MATH_CONVEXHULL_HPP_*/
