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
 
#include "minimizers/MinimizerFactory.hpp"
#include "minimizers/MinimizerCreator.hpp"

namespace oopse {

//initialize instance of MinimizerFactory
MinimizerFactory* MinimizerFactory::instance_ = NULL;

MinimizerFactory::~MinimizerFactory() {
    CreatorMapType::iterator i;
    for (i = creatorMap_.begin(); i != creatorMap_.end(); ++i) {
        delete i->second;
    }
    creatorMap_.clear();
}

bool MinimizerFactory::registerMinimizer(MinimizerCreator* creator) {
    return creatorMap_.insert(
        CreatorMapType::value_type(creator->getIdent(), creator)).second;
}

bool MinimizerFactory::unregisterMinimizer(const std::string& id) {
    return creatorMap_.erase(id) == 1;
}

Minimizer* MinimizerFactory::createMinimizer(const std::string& id, SimInfo* info) {
    CreatorMapType::iterator i = creatorMap_.find(id);
    if (i != creatorMap_.end()) {
        //invoke functor to create object
        return (i->second)->create(info);
    } else {
        return NULL;
    }
}

std::vector<std::string> MinimizerFactory::getIdents() {
    IdentVectorType idents;
    CreatorMapType::iterator i;

    for (i = creatorMap_.begin(); i != creatorMap_.end(); ++i) {
        idents.push_back(i->first);
    }
    
    return idents;
}

std::ostream& operator <<(std::ostream& o, MinimizerFactory& factory) {
    MinimizerFactory::IdentVectorType idents;
    MinimizerFactory::IdentVectorIterator i;

    idents = factory.getIdents();

    o << "Avaliable type identifiers in this factory: " << std::endl;
    for (i = idents.begin(); i != idents.end(); ++i) {
        o << *i << std::endl;
    }

    return o;
}

}

