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
 *  NanoVolume.hpp
 *
 *  Purpose: To calculate convexhull, hull volume and radius
 *  using the CGAL library.
 *
 *  Created by Charles F. Vardeman II on 14 Dec 2006.
 *  @author  Charles F. Vardeman II
 *  @version $Id: NanoVolume.hpp,v 1.5 2008-09-14 01:32:23 chuckv Exp $
 *
 */
#ifndef APPLICATIONS_STATICPROPS_NANOVOLUME_HPP_
#define APPLICATIONS_STATICPROPS_NANOVOLUME_HPP_
#include <vector>
#include "config.h"
#include "math/Vector3.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "applications/staticProps/StaticAnalyser.hpp"

#if defined(HAVE_QHULL) || defined(HAVE_CGAL)
#ifdef HAVE_QHULL
#include "math/ConvexHull.hpp"
#endif

#ifdef HAVE_CGAL
#include "math/AlphaShape.hpp"
#endif
#endif

namespace oopse {
  class NanoVolume : public StaticAnalyser {
  public:
    NanoVolume(SimInfo* info, const std::string& filename, const std::string& sele);
    virtual void process();
    
  private:    
    Snapshot* currentSnapshot_;
    std::string selectionScript_;
    SelectionManager seleMan_;
    SelectionEvaluator evaluator_;
    std::vector<Atom*> theAtoms_;
    int frameCounter_;
    RealType totalVolume_;
    
  };
}
#endif /*NANOVOLUME_HPP_*/
