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
 
/**
 * @file EAM_FF.hpp
 * @author tlin
 * @date 11/15/2004
 * @time 13:43am
 * @version 1.0
 */
#ifndef USETHEFORCE_EAM_FF_HPP
#define USETHEFORCE_EAM_FF_HPP

#include "UseTheForce/ForceField.hpp"
#include "io/SectionParserManager.hpp"
#include "UseTheForce/ForceFieldCreator.hpp"
#include "io/ForceFieldOptions.hpp"
namespace oopse {

  /**
   * @class EAM_FF EAM_FF.hpp "UseTheForce/EAM_FF.hpp"
   * @brief Embedded Atom Method Force Field Parser
   */
  class EAM_FF : public ForceField {
  public:
    EAM_FF();
    virtual void parse(const std::string& filename);

    virtual RealType getRcutFromAtomType(AtomType* at);

    ~EAM_FF();
  private:
    SectionParserManager spMan_;

  };

} //end namespace oopse
#endif //USETHEFORCE_EAM_FF_HPP

