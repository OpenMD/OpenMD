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
 *  sphericalLattice.hpp
 *
 *  Created by Charles F. Vardeman II on 17 Feb 2006.
 *  @author  Charles F. Vardeman II
 *  @version $Id: shapedLattice.hpp,v 1.6 2006-10-17 15:24:29 gezelter Exp $
 *
 */

#ifndef LATTICE_SHAPEDLATTICE_HPP
#define LATTICE_SHAPEDLATTICE_HPP 

#include "math/Vector3.hpp"
#include "lattice/LatticeFactory.hpp"
#include "lattice/Lattice.hpp"
#include "brains/Register.hpp"

namespace oopse{
  
  /**
   * Returns a vector of vector3 position on a lattice truncated 
   * 
   */
   
  class shapedLattice{
  public:
    shapedLattice(RealType latticeConstant, std::string latticeType);
    virtual ~shapedLattice(){};
    /**
     * setGridDimension:  
     * 
     */
    void setGridDimension(Vector3d dimension);
    void setOrigin(Vector3d origin);
    virtual bool isInterior(Vector3d point) =0;
    std::vector<Vector3d> getSites();
    std::vector<Vector3d> getOrientations();
  protected:
    void findSites();
    Vector3d dimension_;
    Vector3d origin_;  
  private:
    bool sitesComputed_;
    std::vector<Vector3d> sites_;
    std::vector<Vector3d> orientations_;
    Lattice *simpleLattice_;
    RealType latticeConstant_;
    std::string latticeType_;
    int beginNx_;
    int beginNy_;
    int beginNz_;
    int endNx_;
    int endNy_;
    int endNz_;
    
  };
}
#endif /* LATTICE_SHAPEDLATTICE_HPP */
