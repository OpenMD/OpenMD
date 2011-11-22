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
 
/**
 * @file ForceField.hpp
 * @author tlin
 * @date 11/04/2004
 * @time 22:51am
 * @version 1.0
 */
  
#ifndef USETHEFORCE_FORCEFIELD_HPP
#define USETHEFORCE_FORCEFIELD_HPP

#include "config.h"
#include <string>
#include <utility>
#include <vector>

#include "io/ifstrstream.hpp"
#include "io/ForceFieldOptions.hpp"
#include "utils/TypeContainer.hpp"
#include "utils/Tuple.hpp"
#include "types/AtomType.hpp"
#include "types/BondType.hpp"
#include "types/BendType.hpp"
#include "types/TorsionType.hpp"
#include "types/InversionType.hpp"
#include "types/NonBondedInteractionType.hpp"
namespace OpenMD {

  /**
   * @class ForceField ForceField.hpp ''UseTheForce/ForceField.hpp"
   * @brief
   */
  class ForceField{

  public:

    typedef TypeContainer<AtomType, 1> AtomTypeContainer;
    typedef TypeContainer<BondType, 2> BondTypeContainer;
    typedef TypeContainer<BendType, 3> BendTypeContainer;
    typedef TypeContainer<TorsionType, 4> TorsionTypeContainer;
    typedef TypeContainer<InversionType, 4> InversionTypeContainer;
    typedef TypeContainer<NonBondedInteractionType, 2> NonBondedInteractionTypeContainer;
    
    ForceField();
    
    virtual ~ForceField() {};

    std::string getForceFieldFileName() {
      return forceFieldFileName_;
    }

    void setForceFieldFileName(const std::string& filename) {
      forceFieldFileName_ = filename;
    }
        
    virtual void parse(const std::string& filename) = 0;  

    AtomType* getAtomType(const std::string &at);
    AtomType* getAtomType(int ident);
    BondType* getBondType(const std::string &at1, const std::string &at2);
    BendType* getBendType(const std::string &at1, const std::string &at2,
                          const std::string &at3);
    TorsionType* getTorsionType(const std::string &at1, const std::string &at2,
                                const std::string &at3, const std::string &at4);
    InversionType* getInversionType(const std::string &at1, 
				    const std::string &at2,
				    const std::string &at3, 
				    const std::string &at4);
    NonBondedInteractionType* getNonBondedInteractionType(const std::string &at1, const std::string &at2);
    
    BondType* getExactBondType(const std::string &at1, const std::string &at2);
    BendType* getExactBendType(const std::string &at1, const std::string &at2,
                               const std::string &at3);
    TorsionType* getExactTorsionType(const std::string &at1, 
                                     const std::string &at2,
                                     const std::string &at3, 
                                     const std::string &at4);
    InversionType* getExactInversionType(const std::string &at1, 
					 const std::string &at2,
					 const std::string &at3, 
					 const std::string &at4);
    NonBondedInteractionType* getExactNonBondedInteractionType(const std::string &at1, const std::string &at2);
    
    
    //avoid make virtual function public
    //Herb Sutter and Andrei Alexandrescu, C++ coding Standards, Addision-Wesley
    virtual RealType getRcutFromAtomType(AtomType* at);
    
    std::string getWildCard() {
      return wildCardAtomTypeName_;
    }
    
    void setWildCard(const std::string& wildCard) {
      wildCardAtomTypeName_ = wildCard;
    }
    
    unsigned int getNAtomType() {
      return atomTypeCont_.size();
    }

    AtomTypeContainer* getAtomTypes() {
      return &atomTypeCont_;
    }
    
    NonBondedInteractionTypeContainer* getNonBondedInteractionTypes() {
      return &nonBondedInteractionTypeCont_;
    }
        
    bool addAtomType(const std::string &at, AtomType* atomType);

    bool replaceAtomType(const std::string &at, AtomType* atomType);

    bool addBondType(const std::string &at1, const std::string &at2, 
                     BondType* bondType);

    bool addBendType(const std::string &at1, const std::string &at2,
		     const std::string &at3, BendType* bendType);

    bool addTorsionType(const std::string &at1, const std::string &at2,
			const std::string &at3, const std::string &at4, 
			TorsionType* torsionType);

    bool addInversionType(const std::string &at1, const std::string &at2,
			  const std::string &at3, const std::string &at4,
			  InversionType* inversionType);
    
    bool addNonBondedInteractionType(const std::string &at1, 
				     const std::string &at2, 
				     NonBondedInteractionType* nbiType);
    
    ifstrstream* openForceFieldFile(const std::string& filename);
    
    ForceFieldOptions& getForceFieldOptions() {return forceFieldOptions_;}

  protected:

    AtomTypeContainer atomTypeCont_;    
    BondTypeContainer bondTypeCont_;
    BendTypeContainer bendTypeCont_;
    TorsionTypeContainer torsionTypeCont_;
    InversionTypeContainer inversionTypeCont_;
    NonBondedInteractionTypeContainer nonBondedInteractionTypeCont_;
    ForceFieldOptions forceFieldOptions_;
    std::map<int, std::string> atypeIdentToName;
    
  private:  
    std::string ffPath_;    
    std::string wildCardAtomTypeName_;    
    std::string forceFieldFileName_;    
    
  };
}//end namespace OpenMD
#endif

