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
 
#ifndef IO_SECTIONPARSERMANAGER_HPP
#define IO_SECTIONPARSERMANAGER_HPP

#include <list>

#include "io/SectionParser.hpp"

namespace OpenMD {

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
