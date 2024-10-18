/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
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
    std::streampos offset;
    bool isActive;
  };

  class SameSectionParserFunctor {
  public:
    SameSectionParserFunctor(const std::string& section) : section_(section) {}

    bool operator()(SectionParserContext context) {
      return context.sectionParser->getSectionName() == section_;
    }

  private:
    std::string section_;
  };
  /**
   * @class SectionParserManager SectionParserManager.hpp
   * "io/SectionParserManager.hpp" SectionParserManager maintains a priority
   * list
   */
  class SectionParserManager {
  public:
    using SectionParserContextList = std::list<SectionParserContext>;
    using iterator                 = SectionParserContextList::iterator;
    using const_iterator           = SectionParserContextList::const_iterator;

    SectionParserManager() : beginPriority_(0), priorityDifference_(100) {}
    ~SectionParserManager();

    void parse(std::istream& input, ForceField& ff);

    void push_front(SectionParser* sp);

    void push_back(SectionParser* sp);

    void insert(SectionParser* sp, int priority);

    const_iterator begin() const { return sectionParsers_.begin(); }

    const_iterator end() const { return sectionParsers_.end(); }

  private:
    iterator findSectionParser(const std::string& sectionName);
    const int beginPriority_;
    int priorityDifference_;

    SectionParserContextList sectionParsers_;
  };
}  // namespace OpenMD

#endif  // IO_SECTIONPARSERMANAGER_HPP
