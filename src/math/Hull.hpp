/* Copyright (c) 2008 The University of Notre Dame. All Rights Reserved.
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
 *  Hull.hpp
 *
 *  Purpose: Provide common interface for computing hulls in oopse.
 *
 *  Created by Charles F. Vardeman II on 27 July 2008.
 *  @author  Charles F. Vardeman II
 *  @version $Id: Hull.hpp,v 1.3 2008-10-21 16:44:00 chuckv Exp $
 *
 */

#ifndef MATH_HULL_HPP_
#define MATH_HULL_HPP_

#include "config.h"
#include "math/Vector3.hpp"
#include "math/Triangle.hpp"
#include "primitives/StuntDouble.hpp"

#include <cassert>
#include <vector>
#include <string>



namespace oopse {
  class Hull {
  public:
    virtual ~Hull(){};
    virtual void computeHull(std::vector<StuntDouble*> bodydoubles)=0;
    virtual RealType getArea()=0; //Total area of Hull
    virtual int getNs()=0;  //Number of Surface Atoms
    virtual int getNMeshElements()=0; //Number of polygons in surface mesh
    virtual RealType getVolume()=0; //Total Volume inclosed by Hull
    virtual std::vector< StuntDouble* > getSurfaceAtoms()=0; //Returns a list of surface atoms
    virtual std::vector<Triangle > getMesh()=0;
    virtual void printHull(const std::string& geomFileName)=0;
  };
}

#endif /*MATH_HULL_HPP_*/
