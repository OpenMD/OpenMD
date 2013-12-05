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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
#include "selection/NameFinder.hpp"
#include "utils/wildcards.hpp"
#include "utils/StringTokenizer.hpp"
#include "primitives/Molecule.hpp"
#include "utils/StringUtils.hpp"
namespace OpenMD {

  TreeNode::~TreeNode(){
    std::map<std::string, TreeNode*>::iterator i;
    for ( i = children.begin(); i != children.end(); ++i) {
      i->second->~TreeNode();
    }
    children.clear();
  }

  NameFinder::NameFinder(SimInfo* info) : info_(info), root_(NULL){
    nObjects_.push_back(info_->getNGlobalAtoms()+info_->getNGlobalRigidBodies());
    nObjects_.push_back(info_->getNGlobalBonds());
    nObjects_.push_back(info_->getNGlobalBends());
    nObjects_.push_back(info_->getNGlobalTorsions());
    nObjects_.push_back(info_->getNGlobalInversions());
    loadNames();
  }

  NameFinder::~NameFinder(){
    delete root_;
  }

  void NameFinder::loadNames() {
    SimInfo::MoleculeIterator mi;
    Molecule::AtomIterator ai;
    Molecule::RigidBodyIterator rbIter;
    Molecule::BondIterator bondIter;
    Molecule::BendIterator bendIter;
    Molecule::TorsionIterator torsionIter;
    Molecule::InversionIterator inversionIter;

    Molecule* mol;
    Atom* atom;
    RigidBody* rb;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Inversion* inversion;    

    root_ = new TreeNode;
    root_->bs.resize(nObjects_);
    root_->bs.setAll(); //
    
    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) {
           
      std::string molName = mol->getMoleculeName();
      TreeNode* molNode = createNode(root_, molName);
        
      for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
	std::string atomName = atom->getType();
	TreeNode* atomNode = createNode(molNode, atomName);
            
	molNode->bs.bitsets_[STUNTDOUBLE].setBitOn(atom->getGlobalIndex());
	atomNode->bs.bitsets_[STUNTDOUBLE].setBitOn(atom->getGlobalIndex());
      }

      for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
           rb = mol->nextRigidBody(rbIter)) {
	std::string rbName = rb->getType();
	TreeNode* rbNode = createNode(molNode, rbName);
            
	molNode->bs.bitsets_[STUNTDOUBLE].setBitOn(rb->getGlobalIndex());
	rbNode->bs.bitsets_[STUNTDOUBLE].setBitOn(rb->getGlobalIndex());

	//create nodes for atoms belong to this rigidbody
	for(atom = rb->beginAtom(ai); atom != NULL; atom = rb->nextAtom(ai)) {
	  std::string rbAtomName = atom->getType();
	  TreeNode* rbAtomNode = createNode(rbNode, rbAtomName);

	  rbAtomNode->bs.bitsets_[STUNTDOUBLE].setBitOn(atom->getGlobalIndex());
	}
      }

      for (bond = mol->beginBond(bondIter); bond != NULL; 
           bond = mol->nextBond(bondIter)) {

	std::string bondName = bond->getName();
	TreeNode* bondNode = createNode(molNode, bondName);

	molNode->bs.bitsets_[BOND].setBitOn(bond->getGlobalIndex());
        bondNode->bs.bitsets_[BOND].setBitOn(bond->getGlobalIndex());

        std::vector<Atom*> atoms = bond->getAtoms();
        std::vector<Atom*>::iterator ai;
        
        for (ai = atoms.begin(); ai != atoms.end(); ++ai) {
          std::string atomName = (*ai)->getType();
          TreeNode* atomNode = createNode(bondNode, atomName);
          atomNode->bs.bitsets_[STUNTDOUBLE].setBitOn((*ai)->getGlobalIndex());
        }
      }
      for (bend = mol->beginBend(bendIter); bend != NULL; 
           bend = mol->nextBend(bendIter)) {

	std::string bendName = bend->getName();
	TreeNode* bendNode = createNode(molNode, bendName);

	molNode->bs.bitsets_[BEND].setBitOn(bend->getGlobalIndex());
        bendNode->bs.bitsets_[BEND].setBitOn(bend->getGlobalIndex());

        std::vector<Atom*> atoms = bend->getAtoms();
        std::vector<Atom*>::iterator ai;
        
        for (ai = atoms.begin(); ai != atoms.end(); ++ai) {
          std::string atomName = (*ai)->getType();
          TreeNode* atomNode = createNode(bendNode, atomName);
          atomNode->bs.bitsets_[STUNTDOUBLE].setBitOn((*ai)->getGlobalIndex());
        }

      }
      for (torsion = mol->beginTorsion(torsionIter); torsion != NULL; 
           torsion = mol->nextTorsion(torsionIter)) {

	std::string torsionName = torsion->getName();
	TreeNode* torsionNode = createNode(molNode, torsionName);

	molNode->bs.bitsets_[TORSION].setBitOn(torsion->getGlobalIndex());
        torsionNode->bs.bitsets_[TORSION].setBitOn(torsion->getGlobalIndex());

        std::vector<Atom*> atoms = torsion->getAtoms();
        std::vector<Atom*>::iterator ai;
        
        for (ai = atoms.begin(); ai != atoms.end(); ++ai) {
          std::string atomName = (*ai)->getType();
          TreeNode* atomNode = createNode(torsionNode, atomName);
          atomNode->bs.bitsets_[STUNTDOUBLE].setBitOn((*ai)->getGlobalIndex());
        }

      }
      for (inversion = mol->beginInversion(inversionIter); inversion != NULL; 
           inversion = mol->nextInversion(inversionIter)) {

	std::string inversionName = inversion->getName();
	TreeNode* inversionNode = createNode(molNode, inversionName);

	molNode->bs.bitsets_[INVERSION].setBitOn(inversion->getGlobalIndex());
        inversionNode->bs.bitsets_[INVERSION].setBitOn(inversion->getGlobalIndex());
        std::vector<Atom*> atoms = inversion->getAtoms();
        std::vector<Atom*>::iterator ai;
        
        for (ai = atoms.begin(); ai != atoms.end(); ++ai) {
          std::string atomName = (*ai)->getType();
          TreeNode* atomNode = createNode(inversionNode, atomName);
          atomNode->bs.bitsets_[STUNTDOUBLE].setBitOn((*ai)->getGlobalIndex());
        }
      }
    }
  }

  TreeNode* NameFinder::createNode(TreeNode* parent, const std::string& name) {
    TreeNode* node;    
    std::map<std::string, TreeNode*>::iterator foundIter;
    foundIter = parent->children.find(name);
    if ( foundIter  == parent->children.end()) {
      node = new TreeNode;
      node->name = name;
      node->bs.resize(nObjects_);
      parent->children.insert(std::make_pair(name, node));
    }else {
      node = foundIter->second;
    }
    return node;
  }

  SelectionSet NameFinder::match(const std::string& name){
    SelectionSet bs(nObjects_);
  
    StringTokenizer tokenizer(name, ".");

    std::vector<std::string> names;
    while(tokenizer.hasMoreTokens()) {
      names.push_back(tokenizer.nextToken());
    }

    int size = names.size();

    switch(size) {
    case 1 :
      //could be molecule name, atom name and rigidbody name
      matchMolecule(names[0], bs);
      matchStuntDouble("*", names[0], bs);
      matchBond("*", names[0], bs);
      matchBend("*", names[0], bs);
      matchTorsion("*", names[0], bs);
      matchInversion("*", names[0], bs);
            
      break;
    case 2:
      //could be molecule.*(include atoms and rigidbodies) or rigidbody.*(atoms belong to rigidbody)

      if (!isInteger(names[1])){
	matchRigidAtoms("*", names[0], names[1], bs);
	matchStuntDouble(names[0], names[1], bs);
      } else {
	int internalIndex = lexi_cast<int>(names[1]);
	if (internalIndex < 0) {
	  std::cerr << names[0] << ". " << names[1] << " is an invalid name" << std::endl;           
	} else {
	  matchInternalIndex(names[0], internalIndex, bs);
	}
      }
            
      break;
    case 3:
      //must be molecule.rigidbody.*
      matchRigidAtoms(names[0], names[1], names[2], bs);
      break;
    default:      
      std::cerr << "invalid name: " << name << std::endl;
      break;           
    }

    return bs; 
  }

  void NameFinder::matchMolecule(const std::string& molName, SelectionSet& bs) {
    std::vector<TreeNode*> molNodes = getMatchedChildren(root_, molName);            
    std::vector<TreeNode*>::iterator i;
    for( i = molNodes.begin(); i != molNodes.end(); ++i ) {
      bs |= (*i)->bs;
    }    
  }

  void NameFinder::matchStuntDouble(const std::string& molName, const std::string& sdName, SelectionSet& bs){
    std::vector<TreeNode*> molNodes = getMatchedChildren(root_, molName);            
    std::vector<TreeNode*>::iterator i;
    for( i = molNodes.begin(); i != molNodes.end(); ++i ) {
      std::vector<TreeNode*> sdNodes = getMatchedChildren(*i, sdName);   
      std::vector<TreeNode*>::iterator j;
      for (j = sdNodes.begin(); j != sdNodes.end(); ++j) {
	bs |= (*j)->bs;
      }
    }

  }

  void NameFinder::matchBond(const std::string& molName, 
                             const std::string& bondName, SelectionSet& bs){
    std::vector<TreeNode*> molNodes = getMatchedChildren(root_, molName);            
    std::vector<TreeNode*>::iterator i;
    for( i = molNodes.begin(); i != molNodes.end(); ++i ) {
      std::vector<TreeNode*> bondNodes = getMatchedChildren(*i, bondName);   
      std::vector<TreeNode*>::iterator j;
      for (j = bondNodes.begin(); j != bondNodes.end(); ++j) {
        bs |= (*j)->bs;
	std::vector<TreeNode*> bondAtomNodes = getAllChildren(*j);
	std::vector<TreeNode*>::iterator k;
	for(k = bondAtomNodes.begin(); k != bondAtomNodes.end(); ++k){
	  bs |= (*k)->bs;
	}
      }
    }
  }

  void NameFinder::matchBend(const std::string& molName, const std::string& bendName,  SelectionSet& bs){
    std::vector<TreeNode*> molNodes = getMatchedChildren(root_, molName);            
    std::vector<TreeNode*>::iterator i;
    for( i = molNodes.begin(); i != molNodes.end(); ++i ) {
      std::vector<TreeNode*> bendNodes = getMatchedChildren(*i, bendName);   
      std::vector<TreeNode*>::iterator j;
      for (j = bendNodes.begin(); j != bendNodes.end(); ++j) {
	std::vector<TreeNode*> bendAtomNodes = getAllChildren(*j);
	std::vector<TreeNode*>::iterator k;
	for(k = bendAtomNodes.begin(); k != bendAtomNodes.end(); ++k){
	  bs |= (*k)->bs;
	}
      }
    }
  }
  void NameFinder::matchTorsion(const std::string& molName, const std::string& torsionName, SelectionSet& bs){
    std::vector<TreeNode*> molNodes = getMatchedChildren(root_, molName);            
    std::vector<TreeNode*>::iterator i;
    for( i = molNodes.begin(); i != molNodes.end(); ++i ) {
      std::vector<TreeNode*> torsionNodes = getMatchedChildren(*i, torsionName);   
      std::vector<TreeNode*>::iterator j;
      for (j = torsionNodes.begin(); j != torsionNodes.end(); ++j) {
	std::vector<TreeNode*> torsionAtomNodes = getAllChildren(*j);
	std::vector<TreeNode*>::iterator k;
	for(k = torsionAtomNodes.begin(); k != torsionAtomNodes.end(); ++k){
	  bs |= (*k)->bs;
	}
      }
    }
  }
  void NameFinder::matchInversion(const std::string& molName, const std::string& inversionName, SelectionSet& bs){
    std::vector<TreeNode*> molNodes = getMatchedChildren(root_, molName);            
    std::vector<TreeNode*>::iterator i;
    for( i = molNodes.begin(); i != molNodes.end(); ++i ) {
      std::vector<TreeNode*> inversionNodes = getMatchedChildren(*i, inversionName);   
      std::vector<TreeNode*>::iterator j;
      for (j = inversionNodes.begin(); j != inversionNodes.end(); ++j) {
	std::vector<TreeNode*> inversionAtomNodes = getAllChildren(*j);
	std::vector<TreeNode*>::iterator k;
	for(k = inversionAtomNodes.begin(); k != inversionAtomNodes.end(); ++k){
	  bs |= (*k)->bs;
	}
      }
    }
  }

  void NameFinder::matchRigidAtoms(const std::string& molName, const std::string& rbName, const std::string& rbAtomName, SelectionSet& bs){
    std::vector<TreeNode*> molNodes = getMatchedChildren(root_, molName);            
    std::vector<TreeNode*>::iterator i;
    for( i = molNodes.begin(); i != molNodes.end(); ++i ) {
      std::vector<TreeNode*> rbNodes = getMatchedChildren(*i, rbName);   
      std::vector<TreeNode*>::iterator j;
      for (j = rbNodes.begin(); j != rbNodes.end(); ++j) {
	std::vector<TreeNode*> rbAtomNodes = getMatchedChildren(*j, rbAtomName);
	std::vector<TreeNode*>::iterator k;
	for(k = rbAtomNodes.begin(); k != rbAtomNodes.end(); ++k){
	  bs |= (*k)->bs;
	}
      }
    }

  }

  std::vector<TreeNode*> NameFinder::getAllChildren(TreeNode* node)  {
    std::vector<TreeNode*> childNodes;
    std::map<std::string, TreeNode*>::iterator i;
    for (i = node->children.begin(); i != node->children.end(); ++i) {
      childNodes.push_back(i->second);
    }
    return childNodes;
  }

  std::vector<TreeNode*> NameFinder::getMatchedChildren(TreeNode* node, const std::string& name) {
    std::vector<TreeNode*> matchedNodes;
    std::map<std::string, TreeNode*>::iterator i;
    for (i = node->children.begin(); i != node->children.end(); ++i) {
      if (isMatched( i->first, name)) {
	matchedNodes.push_back(i->second);
      }
    }

    return matchedNodes;
  }

  bool NameFinder::isMatched(const std::string& str, const std::string& wildcard) {
    return Wildcard::wildcardfit(wildcard.c_str(), str.c_str()) > 0 ? true : false;
  }


  void NameFinder::matchInternalIndex(const std::string& name, int internalIndex, SelectionSet& bs){

    SimInfo::MoleculeIterator mi;
    Molecule* mol;

    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) {
           
      if (isMatched(mol->getMoleculeName(), name) ) {
	int natoms = mol->getNAtoms();
	int nrigidbodies = mol->getNRigidBodies();
	if (internalIndex >= natoms + nrigidbodies) {
	  continue;
	} else if (internalIndex < natoms) {
	  bs.bitsets_[STUNTDOUBLE].setBitOn(mol->getAtomAt(internalIndex)->getGlobalIndex());
	  continue;
	} else if ( internalIndex < natoms + nrigidbodies) {
	  bs.bitsets_[STUNTDOUBLE].setBitOn(mol->getRigidBodyAt(internalIndex - natoms)->getGlobalIndex());
	}
      }
    }
  }

  bool NameFinder::isInteger(const std::string &str) {
    for(unsigned int i = 0; i < str.size(); ++i){
      if (!std::isdigit(str[i])) {
	return false;
      }
    }
    return true;
  }
}
