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
 

#ifndef CONSTRAINTS_CONTRAINTPAIR_HPP
#define CONSTRAINTS_CONTRAINTPAIR_HPP

#include "primitives/StuntDouble.hpp"
#include "constraints/ConstraintElem.hpp"
namespace oopse {


  /**
   * @class ConstraintPair
   * @todo document
   */
  class ConstraintPair {

  public:

    ConstraintPair(ConstraintElem* elem1, ConstraintElem* elem2, RealType len) 
      : consElem1_(elem1), consElem2_(elem2), dist2(len*len) { }

    ~ConstraintPair() {
      delete consElem1_;
      delete consElem2_;
    }
    /** Return the first constraint elemet */
    ConstraintElem* getConsElem1() {return consElem1_;}

    /** Retunr the second constraint element */
    ConstraintElem* getConsElem2() {return consElem2_;}

    bool isMoved() { return consElem1_->getMoved() || consElem2_->getMoved(); }        
    RealType getConsDistSquare() {return dist2;}

  private:
        
    ConstraintElem* consElem1_;
    ConstraintElem* consElem2_;        
    RealType dist2;
  };

}

#endif 
