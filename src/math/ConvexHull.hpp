/* Copyright (c) 2006 The University of Notre Dame. All Rights Reserved.
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
 *  ConvexHull.hpp
 *
 *  Purpose: To calculate convexhull, hull volume and radius
 *  using the CGAL library.
 *
 *  Created by Charles F. Vardeman II on 11 Dec 2006.
 *  @author  Charles F. Vardeman II
 *  @version $Id: ConvexHull.hpp,v 1.3 2007-05-30 18:47:04 chuckv Exp $
 *
 */

#ifndef MATH_CONVEXHULL_HPP_
#define MATH_CONVEXHULL_HPP_

#include "math/Vector3.hpp"

#include <cassert>
#include <vector>
#include <string>
extern "C"
{
#include "QuickHull/qhull.h"
#include "QuickHull/mem.h"
#include "QuickHull/qset.h"
#include "QuickHull/geom.h"
#include "QuickHull/merge.h"
#include "QuickHull/poly.h"
#include "QuickHull/io.h"
#include "QuickHull/stat.h"
}



namespace oopse {
  class ConvexHull {
  public:
    ConvexHull();
    virtual ~ConvexHull() {}
    bool genHull(std::vector<Vector3d> pos);
    //std::vector<Vector3d> getHull();
    RealType getVolume();
    //RealType getRadius();
    // RealType getInscribedRadius();
    void geomviewHull(const std::string& geomFileName);
  protected:
    double volume_;
    double area_;
    int dim_;
    const std::string options_;
    
  };
}

#endif /*MATH_CONVEXHULL_HPP_*/
