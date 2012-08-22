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
 * @file ForceField.cpp
 * @author tlin
 * @date 11/04/2004
 * @time 22:51am
 * @version 1.0
 */
  
#include <algorithm>
#include "brains/ForceField.hpp"
#include "utils/simError.h"

#include "io/OptionSectionParser.hpp"
#include "io/BaseAtomTypesSectionParser.hpp"
#include "io/DirectionalAtomTypesSectionParser.hpp"
#include "io/AtomTypesSectionParser.hpp"
#include "io/BendTypesSectionParser.hpp"
#include "io/BondTypesSectionParser.hpp"
#include "io/ChargeAtomTypesSectionParser.hpp"
#include "io/EAMAtomTypesSectionParser.hpp"
#include "io/FluctuatingChargeAtomTypesSectionParser.hpp"
#include "io/GayBerneAtomTypesSectionParser.hpp"
#include "io/InversionTypesSectionParser.hpp"
#include "io/LennardJonesAtomTypesSectionParser.hpp"
#include "io/MultipoleAtomTypesSectionParser.hpp"
#include "io/NonBondedInteractionsSectionParser.hpp"
#include "io/PolarizableAtomTypesSectionParser.hpp"
#include "io/SCAtomTypesSectionParser.hpp"
#include "io/ShapeAtomTypesSectionParser.hpp"
#include "io/StickyAtomTypesSectionParser.hpp"
#include "io/StickyPowerAtomTypesSectionParser.hpp"
#include "io/TorsionTypesSectionParser.hpp"

#include "types/LennardJonesAdapter.hpp"
#include "types/EAMAdapter.hpp"
#include "types/SuttonChenAdapter.hpp"
#include "types/GayBerneAdapter.hpp"
#include "types/StickyAdapter.hpp"

namespace OpenMD {

  ForceField::ForceField(std::string ffName) { 

    char* tempPath; 
    tempPath = getenv("FORCE_PARAM_PATH");
    
    if (tempPath == NULL) {
      //convert a macro from compiler to a string in c++
      STR_DEFINE(ffPath_, FRC_PATH );
    } else {
      ffPath_ = tempPath;
    }

    setForceFieldFileName(ffName + ".frc");

    /**
     * The order of adding section parsers is important.
     *
     * OptionSectionParser must come first to set options for other
     * parsers
     * 
     * DirectionalAtomTypesSectionParser should be added before
     * AtomTypesSectionParser, and these two section parsers will
     * actually create "real" AtomTypes (AtomTypesSectionParser will
     * create AtomType and DirectionalAtomTypesSectionParser will
     * create DirectionalAtomType, which is a subclass of AtomType and
     * should come first). 
     *
     * Other AtomTypes Section Parsers will not create the "real"
     * AtomType, they only add and set some attributes of the AtomType
     * (via the Adapters). Thus ordering of these is not important.
     * AtomTypesSectionParser should be added before other atom type
     *
     * The order of BondTypesSectionParser, BendTypesSectionParser and
     * TorsionTypesSectionParser, etc. are not important.
     */

    spMan_.push_back(new OptionSectionParser(forceFieldOptions_));
    spMan_.push_back(new BaseAtomTypesSectionParser());
    spMan_.push_back(new DirectionalAtomTypesSectionParser(forceFieldOptions_));
    spMan_.push_back(new AtomTypesSectionParser());

    spMan_.push_back(new LennardJonesAtomTypesSectionParser(forceFieldOptions_));
    spMan_.push_back(new ChargeAtomTypesSectionParser(forceFieldOptions_));
    spMan_.push_back(new MultipoleAtomTypesSectionParser(forceFieldOptions_));
    spMan_.push_back(new FluctuatingChargeAtomTypesSectionParser(forceFieldOptions_));
    spMan_.push_back(new PolarizableAtomTypesSectionParser(forceFieldOptions_));
    spMan_.push_back(new GayBerneAtomTypesSectionParser(forceFieldOptions_));
    spMan_.push_back(new EAMAtomTypesSectionParser(forceFieldOptions_));
    spMan_.push_back(new SCAtomTypesSectionParser(forceFieldOptions_));
    spMan_.push_back(new ShapeAtomTypesSectionParser(forceFieldOptions_));
    spMan_.push_back(new StickyAtomTypesSectionParser(forceFieldOptions_));
    spMan_.push_back(new StickyPowerAtomTypesSectionParser(forceFieldOptions_));

    spMan_.push_back(new BondTypesSectionParser(forceFieldOptions_));
    spMan_.push_back(new BendTypesSectionParser(forceFieldOptions_));
    spMan_.push_back(new TorsionTypesSectionParser(forceFieldOptions_));
    spMan_.push_back(new InversionTypesSectionParser(forceFieldOptions_));

    spMan_.push_back(new NonBondedInteractionsSectionParser(forceFieldOptions_));    
  }

  void ForceField::parse(const std::string& filename) {
    ifstrstream* ffStream;

    ffStream = openForceFieldFile(filename);

    spMan_.parse(*ffStream, *this);

    ForceField::AtomTypeContainer::MapTypeIterator i;
    AtomType* at;

    for (at = atomTypeCont_.beginType(i); at != NULL; 
         at = atomTypeCont_.nextType(i)) {

      // useBase sets the responsibilities, and these have to be done 
      // after the atomTypes and Base types have all been scanned:

      std::vector<AtomType*> ayb = at->allYourBase();      
      if (ayb.size() > 1) {
        for (int j = ayb.size()-1; j > 0; j--) {
          
          ayb[j-1]->useBase(ayb[j]);

        }
      }
    }

    delete ffStream;
  }

  /**
   * getAtomType by string
   *
   * finds the requested atom type in this force field using the string
   * name of the atom type.
   */
  AtomType* ForceField::getAtomType(const std::string &at) {
    std::vector<std::string> keys;
    keys.push_back(at);
    return atomTypeCont_.find(keys);
  }

  /**
   * getAtomType by ident
   *
   * finds the requested atom type in this force field using the
   * integer ident instead of the string name of the atom type.
   */
  AtomType* ForceField::getAtomType(int ident) {   
    std::string at = atypeIdentToName.find(ident)->second;
    return getAtomType(at);
  }

  BondType* ForceField::getBondType(const std::string &at1, 
				    const std::string &at2) {
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    

    //try exact match first
    BondType* bondType = bondTypeCont_.find(keys);
    if (bondType) {
      return bondType;
    } else {
      AtomType* atype1;
      AtomType* atype2;
      std::vector<std::string> at1key;
      at1key.push_back(at1);
      atype1 = atomTypeCont_.find(at1key);
  
      std::vector<std::string> at2key;
      at2key.push_back(at2);
      atype2 = atomTypeCont_.find(at2key);

      // query atom types for their chains of responsibility
      std::vector<AtomType*> at1Chain = atype1->allYourBase();
      std::vector<AtomType*> at2Chain = atype2->allYourBase();

      std::vector<AtomType*>::iterator i;
      std::vector<AtomType*>::iterator j;

      int ii = 0;
      int jj = 0;
      int bondTypeScore;

      std::vector<std::pair<int, std::vector<std::string> > > foundBonds;

      for (i = at1Chain.begin(); i != at1Chain.end(); i++) {
	jj = 0;
	for (j = at2Chain.begin(); j != at2Chain.end(); j++) {

	  bondTypeScore = ii + jj;

	  std::vector<std::string> myKeys;
	  myKeys.push_back((*i)->getName());
	  myKeys.push_back((*j)->getName());

	  BondType* bondType = bondTypeCont_.find(myKeys);
	  if (bondType) {
	    foundBonds.push_back(std::make_pair(bondTypeScore, myKeys));
	  }
	  jj++;
	}
	ii++;
      }


      if (foundBonds.size() > 0) {
        // sort the foundBonds by the score:
        std::sort(foundBonds.begin(), foundBonds.end());
     
        int bestScore = foundBonds[0].first;
        std::vector<std::string> theKeys = foundBonds[0].second;
        
        BondType* bestType = bondTypeCont_.find(theKeys);
        
        return bestType;
      } else {
        //if no exact match found, try wild card match
        return bondTypeCont_.find(keys, wildCardAtomTypeName_);      
      }
    }
  }
  
  BendType* ForceField::getBendType(const std::string &at1, 
				    const std::string &at2,
				    const std::string &at3) {
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    keys.push_back(at3);    

    //try exact match first
    BendType* bendType = bendTypeCont_.find(keys);
    if (bendType) {
      return bendType;
    } else {

      AtomType* atype1;
      AtomType* atype2;
      AtomType* atype3;
      std::vector<std::string> at1key;
      at1key.push_back(at1);
      atype1 = atomTypeCont_.find(at1key);
  
      std::vector<std::string> at2key;
      at2key.push_back(at2);
      atype2 = atomTypeCont_.find(at2key);

      std::vector<std::string> at3key;
      at3key.push_back(at3);
      atype3 = atomTypeCont_.find(at3key);

      // query atom types for their chains of responsibility
      std::vector<AtomType*> at1Chain = atype1->allYourBase();
      std::vector<AtomType*> at2Chain = atype2->allYourBase();
      std::vector<AtomType*> at3Chain = atype3->allYourBase();

      std::vector<AtomType*>::iterator i;
      std::vector<AtomType*>::iterator j;
      std::vector<AtomType*>::iterator k;

      int ii = 0;
      int jj = 0;
      int kk = 0;
      int IKscore;

      std::vector<tuple3<int, int, std::vector<std::string> > > foundBends;

      for (j = at2Chain.begin(); j != at2Chain.end(); j++) {
	ii = 0;
	for (i = at1Chain.begin(); i != at1Chain.end(); i++) {
	  kk = 0;
	  for (k = at3Chain.begin(); k != at3Chain.end(); k++) {
	  
	    IKscore = ii + kk;

	    std::vector<std::string> myKeys;
	    myKeys.push_back((*i)->getName());
	    myKeys.push_back((*j)->getName());
	    myKeys.push_back((*k)->getName());

	    BendType* bendType = bendTypeCont_.find(myKeys);
	    if (bendType) { 
	      foundBends.push_back( make_tuple3(jj, IKscore, myKeys) );
	    }
	    kk++;
	  }
	  ii++;
	}
	jj++;
      }
      
      if (foundBends.size() > 0) {
        std::sort(foundBends.begin(), foundBends.end());
        int jscore = foundBends[0].first;
        int ikscore = foundBends[0].second;
        std::vector<std::string> theKeys = foundBends[0].third;       
        
        BendType* bestType = bendTypeCont_.find(theKeys);  
        return bestType;
      } else {        
	//if no exact match found, try wild card match
	return bendTypeCont_.find(keys, wildCardAtomTypeName_);      
      }
    }
  }

  TorsionType* ForceField::getTorsionType(const std::string &at1, 
					  const std::string &at2,
					  const std::string &at3, 
					  const std::string &at4) {
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    keys.push_back(at3);    
    keys.push_back(at4);    


    //try exact match first
    TorsionType* torsionType = torsionTypeCont_.find(keys);
    if (torsionType) {
      return torsionType;
    } else {

      AtomType* atype1;
      AtomType* atype2;
      AtomType* atype3;
      AtomType* atype4;
      std::vector<std::string> at1key;
      at1key.push_back(at1);
      atype1 = atomTypeCont_.find(at1key);
  
      std::vector<std::string> at2key;
      at2key.push_back(at2);
      atype2 = atomTypeCont_.find(at2key);

      std::vector<std::string> at3key;
      at3key.push_back(at3);
      atype3 = atomTypeCont_.find(at3key);

      std::vector<std::string> at4key;
      at4key.push_back(at4);
      atype4 = atomTypeCont_.find(at4key);

      // query atom types for their chains of responsibility
      std::vector<AtomType*> at1Chain = atype1->allYourBase();
      std::vector<AtomType*> at2Chain = atype2->allYourBase();
      std::vector<AtomType*> at3Chain = atype3->allYourBase();
      std::vector<AtomType*> at4Chain = atype4->allYourBase();

      std::vector<AtomType*>::iterator i;
      std::vector<AtomType*>::iterator j;
      std::vector<AtomType*>::iterator k;
      std::vector<AtomType*>::iterator l;

      int ii = 0;
      int jj = 0;
      int kk = 0;
      int ll = 0;
      int ILscore;
      int JKscore;

      std::vector<tuple3<int, int, std::vector<std::string> > > foundTorsions;

      for (j = at2Chain.begin(); j != at2Chain.end(); j++) {
	kk = 0;
	for (k = at3Chain.begin(); k != at3Chain.end(); k++) {
	  ii = 0;	
	  for (i = at1Chain.begin(); i != at1Chain.end(); i++) {
	    ll = 0;
	    for (l = at4Chain.begin(); l != at4Chain.end(); l++) {
	  
	      ILscore = ii + ll;
	      JKscore = jj + kk;

	      std::vector<std::string> myKeys;
	      myKeys.push_back((*i)->getName());
	      myKeys.push_back((*j)->getName());
	      myKeys.push_back((*k)->getName());
	      myKeys.push_back((*l)->getName());

	      TorsionType* torsionType = torsionTypeCont_.find(myKeys);
	      if (torsionType) { 
		foundTorsions.push_back( make_tuple3(JKscore, ILscore, myKeys) );
	      }
	      ll++;
	    }
	    ii++;
	  }
	  kk++;
	}
	jj++;
      }
      
      if (foundTorsions.size() > 0) {
        std::sort(foundTorsions.begin(), foundTorsions.end());
        int jkscore = foundTorsions[0].first;
        int ilscore = foundTorsions[0].second;
        std::vector<std::string> theKeys = foundTorsions[0].third;
        
        TorsionType* bestType = torsionTypeCont_.find(theKeys);
        return bestType;
      } else {
	//if no exact match found, try wild card match
	return torsionTypeCont_.find(keys, wildCardAtomTypeName_);
      }
    }
  }

  InversionType* ForceField::getInversionType(const std::string &at1, 
					      const std::string &at2,
					      const std::string &at3, 
					      const std::string &at4) {
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    keys.push_back(at3);    
    keys.push_back(at4);    

    //try exact match first
    InversionType* inversionType = inversionTypeCont_.permutedFindSkippingFirstElement(keys);
    if (inversionType) {
      return inversionType;
    } else {
      
      AtomType* atype1;
      AtomType* atype2;
      AtomType* atype3;
      AtomType* atype4;
      std::vector<std::string> at1key;
      at1key.push_back(at1);
      atype1 = atomTypeCont_.find(at1key);
      
      std::vector<std::string> at2key;
      at2key.push_back(at2);
      atype2 = atomTypeCont_.find(at2key);
      
      std::vector<std::string> at3key;
      at3key.push_back(at3);
      atype3 = atomTypeCont_.find(at3key);
      
      std::vector<std::string> at4key;
      at4key.push_back(at4);
      atype4 = atomTypeCont_.find(at4key);

      // query atom types for their chains of responsibility
      std::vector<AtomType*> at1Chain = atype1->allYourBase();
      std::vector<AtomType*> at2Chain = atype2->allYourBase();
      std::vector<AtomType*> at3Chain = atype3->allYourBase();
      std::vector<AtomType*> at4Chain = atype4->allYourBase();

      std::vector<AtomType*>::iterator i;
      std::vector<AtomType*>::iterator j;
      std::vector<AtomType*>::iterator k;
      std::vector<AtomType*>::iterator l;

      int ii = 0;
      int jj = 0;
      int kk = 0;
      int ll = 0;
      int Iscore;
      int JKLscore;
      
      std::vector<tuple3<int, int, std::vector<std::string> > > foundInversions;
      
      for (j = at2Chain.begin(); j != at2Chain.end(); j++) {
	kk = 0;
	for (k = at3Chain.begin(); k != at3Chain.end(); k++) {
	  ii = 0;	
	  for (i = at1Chain.begin(); i != at1Chain.end(); i++) {
	    ll = 0;
	    for (l = at4Chain.begin(); l != at4Chain.end(); l++) {
	      
	      Iscore = ii;
	      JKLscore = jj + kk + ll;
	      
	      std::vector<std::string> myKeys;
	      myKeys.push_back((*i)->getName());
	      myKeys.push_back((*j)->getName());
	      myKeys.push_back((*k)->getName());
	      myKeys.push_back((*l)->getName());
	      
	      InversionType* inversionType = inversionTypeCont_.permutedFindSkippingFirstElement(myKeys);
	      if (inversionType) { 
		foundInversions.push_back( make_tuple3(Iscore, JKLscore, myKeys) );
	      }
	      ll++;
	    }
	    ii++;
	  }
	  kk++;
	}
	jj++;
      }
         
      if (foundInversions.size() > 0) {
        std::sort(foundInversions.begin(), foundInversions.end());
        int iscore = foundInversions[0].first;
        int jklscore = foundInversions[0].second;
        std::vector<std::string> theKeys = foundInversions[0].third;
        
        InversionType* bestType = inversionTypeCont_.permutedFindSkippingFirstElement(theKeys);
        return bestType;
      } else {
	//if no exact match found, try wild card match
	return inversionTypeCont_.find(keys, wildCardAtomTypeName_);
      }
    }
  }
  
  NonBondedInteractionType* ForceField::getNonBondedInteractionType(const std::string &at1, const std::string &at2) {
    
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    
    //try exact match first
    NonBondedInteractionType* nbiType = nonBondedInteractionTypeCont_.find(keys);
    if (nbiType) {
      return nbiType;
    } else {
      AtomType* atype1;
      AtomType* atype2;
      std::vector<std::string> at1key;
      at1key.push_back(at1);
      atype1 = atomTypeCont_.find(at1key);
      
      std::vector<std::string> at2key;
      at2key.push_back(at2);
      atype2 = atomTypeCont_.find(at2key);
      
      // query atom types for their chains of responsibility
      std::vector<AtomType*> at1Chain = atype1->allYourBase();
      std::vector<AtomType*> at2Chain = atype2->allYourBase();
      
      std::vector<AtomType*>::iterator i;
      std::vector<AtomType*>::iterator j;
      
      int ii = 0;
      int jj = 0;
      int nbiTypeScore;
      
      std::vector<std::pair<int, std::vector<std::string> > > foundNBI;
      
      for (i = at1Chain.begin(); i != at1Chain.end(); i++) {
        jj = 0;
        for (j = at2Chain.begin(); j != at2Chain.end(); j++) {
          
          nbiTypeScore = ii + jj;
          
          std::vector<std::string> myKeys;
          myKeys.push_back((*i)->getName());
          myKeys.push_back((*j)->getName());
          
          NonBondedInteractionType* nbiType = nonBondedInteractionTypeCont_.find(myKeys);
          if (nbiType) {
            foundNBI.push_back(std::make_pair(nbiTypeScore, myKeys));
          }
          jj++;
        }
        ii++;
      }
      
      
      if (foundNBI.size() > 0) {
        // sort the foundNBI by the score:
        std::sort(foundNBI.begin(), foundNBI.end());
        
        int bestScore = foundNBI[0].first;
        std::vector<std::string> theKeys = foundNBI[0].second;
        
        NonBondedInteractionType* bestType = nonBondedInteractionTypeCont_.find(theKeys);        
        return bestType;
      } else {
        //if no exact match found, try wild card match
        return nonBondedInteractionTypeCont_.find(keys, wildCardAtomTypeName_);
      }
    }
  }
  
  BondType* ForceField::getExactBondType(const std::string &at1, 
					 const std::string &at2){ 
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    return bondTypeCont_.find(keys);
  }
  
  BendType* ForceField::getExactBendType(const std::string &at1, 
					 const std::string &at2,
					 const std::string &at3){ 
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    keys.push_back(at3);    
    return bendTypeCont_.find(keys);
  }
  
  TorsionType* ForceField::getExactTorsionType(const std::string &at1, 
					       const std::string &at2,
					       const std::string &at3, 
					       const std::string &at4){ 
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    keys.push_back(at3);    
    keys.push_back(at4);   
    return torsionTypeCont_.find(keys);
  }
  
  InversionType* ForceField::getExactInversionType(const std::string &at1, 
						   const std::string &at2,
						   const std::string &at3, 
						   const std::string &at4){ 
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    keys.push_back(at3);    
    keys.push_back(at4);   
    return inversionTypeCont_.find(keys);
  }
  
  NonBondedInteractionType* ForceField::getExactNonBondedInteractionType(const std::string &at1, const std::string &at2){ 
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    return nonBondedInteractionTypeCont_.find(keys);
  }
  

  bool ForceField::addAtomType(const std::string &at, AtomType* atomType) {
    std::vector<std::string> keys;
    keys.push_back(at);
    atypeIdentToName[atomType->getIdent()] = at;
    return atomTypeCont_.add(keys, atomType);
  }

  bool ForceField::replaceAtomType(const std::string &at, AtomType* atomType) {
    std::vector<std::string> keys;
    keys.push_back(at);
    atypeIdentToName[atomType->getIdent()] = at;
    return atomTypeCont_.replace(keys, atomType);
  }

  bool ForceField::addBondType(const std::string &at1, const std::string &at2,
			       BondType* bondType) {
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    return bondTypeCont_.add(keys, bondType);    
  }
  
  bool ForceField::addBendType(const std::string &at1, const std::string &at2,
			       const std::string &at3, BendType* bendType) {
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    keys.push_back(at3);    
    return bendTypeCont_.add(keys, bendType);
  }
  
  bool ForceField::addTorsionType(const std::string &at1, 
				  const std::string &at2,
				  const std::string &at3, 
				  const std::string &at4, 
				  TorsionType* torsionType) {
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    keys.push_back(at3);    
    keys.push_back(at4);    
    return torsionTypeCont_.add(keys, torsionType);
  }

  bool ForceField::addInversionType(const std::string &at1, 
				    const std::string &at2,
				    const std::string &at3, 
				    const std::string &at4, 
				    InversionType* inversionType) {
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    keys.push_back(at3);    
    keys.push_back(at4);    
    return inversionTypeCont_.add(keys, inversionType);
  }
  
  bool ForceField::addNonBondedInteractionType(const std::string &at1, 
					       const std::string &at2, 
					       NonBondedInteractionType* nbiType) {
    std::vector<std::string> keys;
    keys.push_back(at1);
    keys.push_back(at2);    
    return nonBondedInteractionTypeCont_.add(keys, nbiType);
  }
  
  RealType ForceField::getRcutFromAtomType(AtomType* at) {
    RealType rcut(0.0);
    
    LennardJonesAdapter lja = LennardJonesAdapter(at);
    if (lja.isLennardJones()) {
      rcut = 2.5 * lja.getSigma();
    }
    EAMAdapter ea = EAMAdapter(at);
    if (ea.isEAM()) {
      rcut = max(rcut, ea.getRcut());
    }
    SuttonChenAdapter sca = SuttonChenAdapter(at);
    if (sca.isSuttonChen()) {
      rcut = max(rcut, 2.0 * sca.getAlpha());
    }
    GayBerneAdapter gba = GayBerneAdapter(at);
    if (gba.isGayBerne()) {
      rcut = max(rcut, 2.5 * sqrt(2.0) * max(gba.getD(), gba.getL()));
    }
    StickyAdapter sa = StickyAdapter(at);
    if (sa.isSticky()) {
      rcut = max(rcut, max(sa.getRu(), sa.getRup()));
    }

    return rcut;    
  }
  

  ifstrstream* ForceField::openForceFieldFile(const std::string& filename) {
    std::string forceFieldFilename(filename);
    ifstrstream* ffStream = new ifstrstream();
    
    //try to open the force filed file in current directory first    
    ffStream->open(forceFieldFilename.c_str());
    if(!ffStream->is_open()){

      forceFieldFilename = ffPath_ + "/" + forceFieldFilename;
      ffStream->open( forceFieldFilename.c_str() );

      //if current directory does not contain the force field file,
      //try to open it in the path        
      if(!ffStream->is_open()){

	sprintf( painCave.errMsg,
		 "Error opening the force field parameter file:\n"
		 "\t%s\n"
		 "\tHave you tried setting the FORCE_PARAM_PATH environment "
		 "variable?\n",
		 forceFieldFilename.c_str() );
	painCave.severity = OPENMD_ERROR;
	painCave.isFatal = 1;
	simError();
      }
    }  
    return ffStream;
  }

} //end namespace OpenMD
