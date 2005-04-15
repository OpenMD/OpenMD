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
 
#ifndef IO_SECTIONPARSERMANAGER_HPP
#define IO_SECTIONPARSERMANAGER_HPP

#include <list>

#include "io/SectionParser.hpp"

namespace oopse {

  struct SectionParserContext {
    int priority;
    SectionParser* sectionParser;
    int lineNo;
    std::streampos  offset;
    bool isActive;
  };

  class SameSectionParserFunctor {
  public:
    SameSectionParserFunctor(const std::string section) : section_(section) {}

    bool operator()(SectionParserContext context) {
      return context.sectionParser->getSectionName() == section_;
    }
        
  private:
    std::string section_;
  };
  /**
   * @class SectionParserManager SectionParserManager.hpp "io/SectionParserManager.hpp"
   * SectionParserManager maintains a priority list
   */
  class SectionParserManager {

  public:
    typedef std::list<SectionParserContext> SectionParserContextList;
    typedef SectionParserContextList::iterator iterator;
    typedef SectionParserContextList::const_iterator const_iterator;

    SectionParserManager() : beginPriority_(0), priorityDifference_(100) {}
    ~SectionParserManager();

    void parse(std::istream& input, ForceField&  ff);
        
    void push_front(SectionParser* sp);

    void push_back(SectionParser* sp);
        
    void insert(SectionParser* sp, int priority);
        
    const_iterator begin() const {
      return sectionParsers_.begin();
    }

    const_iterator end() const{
      return sectionParsers_.end();
    }
        
  private:
    iterator findSectionParser(const std::string& sectionName);
    const int beginPriority_;
    int priorityDifference_;
        
    SectionParserContextList sectionParsers_;
  };

}
#endif //IO_SECTIONPARSERMANAGER_HPP
