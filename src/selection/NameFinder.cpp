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
#include "selection/NameFinder.hpp"
#include "utils/wildcards.hpp"
namespace oopse {

TreeNode::~TreeNode(){
    std::map<std::string, TreeNode*>::iterator i;
    for ( i = children.begin(); i != children.end(); ++i) {
        i->second->~TreeNode();
    }
    children.clear();
}


NameFinder::NameFinder(SimInfo* info) : info_(info), root_(NULL){
    nStuntDouble_ = info_->getNGlobalAtoms() + info_->getNGlobalRigidBodies()
    loadNames();
}


NameFinder::~NameFinder(){
    delete root_;
}

void NameFinder::loadNames() {

    std::map<std::string, TreeNode*>::iterator foundIter;
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::AtomIterator ai;
    Atom* atom;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;

    root_ = new TreeNode;
    root_->bs.resize(nStuntDobule_);
    root_->bs.setAll(); //
    
    for (mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)) {
        TreeNode* currentMolNode;            
        std::string molName = mol->getMoleculeName();

        foundIter = root_->children.find(molName);
        if ( foundIter  == root_->children.end()) {
            currentMolNode = new TreeNode;
            currentMolNode->name = molName;
            currentMolNode->bs.resize(nStuntDouble_);
        }else {
            currentMolNode = i->second;
        }
        
        for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
            std::string atomName = atom->getType();
            TreeNode* currentAtomNode;
            foundIter = currentMolNode->children.find(molName);
            if (foundIter == currentMolNode->children.end()) {
                currentAtomNode = new TreeNode;
                currentAtomNode->name = atomName;
                currentAtomNode->bs.resize(nStuntDouble_);
            } else {
                currentAtomNode = foundIter->second;
            }
            currentMolNode->bs.setBitOn(atom->getGlobalIndex());
            currentAtomNode->bs.setBitOn(atom->getGlobalIndex());
        }

        for (rb = mol->beginRigidBody(rbIter); rb != NULL; rb = mol->nextRigidBody(rbIter)) {
            std::string rbName = atom->getType();
            TreeNode* currentRbNode;
            foundIter = currentMolNode->children.find(molName);
            if (foundIter == currentMolNode->children.end()) {
                currentRbNode = new TreeNode;
                currentRbNode->name = rbName;
                currentRbNode->bs.resize(nStuntDouble_);
            } else {
                currentRbNode = foundIter->second;
            }
            
            currentMolNode->bs.setBitOn(rb->getGlobalIndex());
            currentRbNode->bs.setBitOn(rb->getGlobalIndex());

            //create nodes for atoms belong to this rigidbody
            for(atom = rb->beginAtom(ai); rb != NULL; atom = rb->nextAtom(ai)) {
                std::string rbAtomName = atom->getType();
                TreeNode* currentRbAtomNode;
                foundIter = currentRbNode->children.find(molName);
                if (foundIter == currentRbNode->children.end()) {
                    currentRbAtomNode = new TreeNode;
                    currentRbAtomNode->name = rbAtomName;
                    currentRbAtomNode->bs.resize(nStuntDouble_);
                } else {
                    currentRbAtomNode = foundIter->second;
                }
                currentRbAtomNode->bs.setBitOn(atom->getGlobalIndex());
            }

        }
        
    }    

    std::map<std::string, TreeNode*>::iterator i;
    for( i = root_->children.begin(); i != ; ++i){
        i->bs = 
    }
}

bool NameFinder::match(const std::string& name, BitSet& bs){

    bool error = true;
    StringTokenizer tokenizer(name, ".");

    std::vector<std::string> names;
    while(tokenizer.hasMoreTokens()) {
        names.push_back(tokenizer.nextToken());
    }

    int size = names.size();
    switch(size) {
        case 1 :
            //could be molecule name, atom name and rigidbody name
            if (names[0] == "*"){
                //if all molecules are selected, we don't need to do the matching, just set all of the bits
                bs.setAll();
            } else{
                matchMolecule(name[0]);
                matchStuntDouble("*", names[0]);
            } 
            
            break;
        case 2:
            //could be molecule.*(include atoms and rigidbodies) or rigidbody.*(atoms belong to rigidbody)
            matchRigidAtoms("*", names[0], names[1], bs);
            matchStuntDouble(names[0], names[1]);
            
            break;
        case 3:
            //must be molecule.rigidbody.*
            matchRigidAtoms(names[0], names[1], names[2], bs)
            break;
        default:            
            break;           
    }

    return matched;
}

void NameFinder::matchMolecule(const std::string& molName, BitSet& bs) {
    std::vector<TreeNode*> molNodes = getMatchedChildren(root_, molName);            
    std::vector<TreeNode*>::iterator i;
    for( i = molNodes.begin(); i != molNodes.end(); ++i ) {
        bs |= i->bs;
    }    
}

void NameFinder::matchStuntDouble(const std::string& molName, const std::string& sdName, BitSet& bs){
    std::vector<TreeNode*> molNodes = getMatchedChildren(root_, molName);            
    std::vector<TreeNode*>::iterator i;
    for( i = molNodes.begin(); i != molNodes.end(); ++i ) {
        std::vector<TreeNode*> sdNodes = getMatchedChildren(*i, sdName);   
        std::vector<TreeNode*>::iterator j;
        for (j = sdNodes.begin(); j != sdNodes.end(); ++j) {
            bs |= j->bs;
        }
    }

}

void NameFinder::matchRigidAtoms(const std::string& molName, const std::string& rbName, const std::string& rbAtomName, BitSet& bs){
    std::vector<TreeNode*> molNodes = getMatchedChildren(root_, molName);            
    std::vector<TreeNode*>::iterator i;
    for( i = molNodes.begin(); i != molNodes.end(); ++i ) {
        std::vector<TreeNode*> rbNodes = getMatchedChildren(*i, rbName);   
        std::vector<TreeNode*>::iterator j;
        for (j = rbNodes.begin(); j != rbNodes.end(); ++j) {
            std::vector<TreeNode*> rbAtomNodes = getMatchedChildren(*j, rbAtomName);
            std::vector<TreeNode*>::iterator k;
            for(k = rbAtomNodes.begin(); k != rbAtomNodes.end(); ++k){
                bs |= k->bs;
            }
        }
    }

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
    return Wildcard::wildcardfit (wildcard.c_str(), str.c_str());
}

}
