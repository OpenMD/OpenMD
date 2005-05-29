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
#include <sstream> 
#include <fstream>
#include "visitors/AtomNameVisitor.hpp"
#include "utils/Trim.hpp"
#include "utils/StringTokenizer.hpp"


namespace oopse {
AtomNameVisitor::AtomNameVisitor(SimInfo* info) : BaseVisitor(), info_(info) {
    visitorName = "AtomNameVisitor";
    std::stringstream ss(defaultAtomTypeTable);
    readAtomTypes(ss);
}

AtomNameVisitor::AtomNameVisitor(SimInfo* info, const std::string& atomTypeFile) : BaseVisitor(), info_(info) {
    visitorName = "AtomNameVisitor";
    std::ifstream ifs(atomTypeFile.c_str());
    readAtomTypes(ifs);
}


void AtomNameVisitor::visitAtom(Atom* atom) {
    AtomData* atomData;
    GenericData* data;    
    bool haveAtomData;
    
    data = atom->getPropertyByName("ATOMDATA");
    if(data != NULL){
      atomData = dynamic_cast<AtomData*>(data);  
      if(atomData == NULL){
	std::cerr << "can not get Atom Data from " << atom->getType() << std::endl;
	atomData = new AtomData; 
	haveAtomData = false;      
      } else {
	haveAtomData = true;
      }
    } else {
      atomData = new AtomData;
      haveAtomData = false;
    }

    std::vector<AtomInfo*>::iterator i;
    for (AtomInfo* atomInfo = atomData->beginAtomInfo(i); atomInfo != NULL; atomInfo = atomData->nextAtomInfo(i)) {
        atomInfo->atomTypeName = getBaseAtomTypeName(atomInfo->atomTypeName);        
    }
    
}

void AtomNameVisitor::visit(RigidBody* rb) {
    std::vector<Atom*>::iterator i;

  for (Atom* atom = rb->beginAtom(i); atom != NULL; atom = rb->nextAtom(i)) {
    visit(atom);
  }
}


void AtomNameVisitor::readAtomTypes(std::istream& is) {
    const int bufferSize = 65535;
    char buffer[bufferSize];
    std::string line;
    while(is.getline(buffer, bufferSize)) {
      line = trimLeftCopy(buffer);
      //a line begins with "//" is comment
      // let's also call lines starting with # and ! as comments
     if ( line.empty() || 
                  (line.size() >= 2 && line[0] == '/' && line[1] == '/') ||
                  (line.size() >= 1 && line[0] == '#') || 
                  (line.size() >= 1 && line[0] == '!') ) {
        continue;
      } else {
        StringTokenizer tokenizer(line);
        if (tokenizer.countTokens() >= 2) {
            std::string atomName = tokenizer.nextToken();
            std::string baseAtomName = tokenizer.nextToken();
            atomNames_.insert(MapType::value_type(atomName, baseAtomName));
        }
      }
    }

}

std::string AtomNameVisitor::getBaseAtomTypeName(const std::string& atomTypeName) {
    MapType::iterator i;
    i = atomNames_.find(atomTypeName);
    return i != atomNames_.end() ? i->second : atomTypeName;
}


const std::string AtomNameVisitor::toString() {
    char   buffer[65535];
    std::string result;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;

    sprintf(buffer, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    sprintf(buffer,
            "Visitor Description: print base atom types\n");
    result += buffer;

    sprintf(buffer,
            "------------------------------------------------------------------\n");
    result += buffer;

    return result;
}

}
