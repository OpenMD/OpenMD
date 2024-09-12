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
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */
#include "io/SectionParserManager.hpp"

#include <algorithm>
#include <cstdio>
#include <stack>

#include "utils/Trim.hpp"
#include "utils/StringUtils.hpp"
#include "utils/simError.h"

namespace OpenMD {

  SectionParserManager::~SectionParserManager() {
    SectionParserManager::iterator i;
    for (i = sectionParsers_.begin(); i != sectionParsers_.end(); ++i) {
      delete (i->sectionParser);
    }
    sectionParsers_.clear();
  }

  void SectionParserManager::parse(std::istream& input, ForceField& ff) {
    // reset active flags
    SectionParserManager::iterator i;
    for (i = sectionParsers_.begin(); i != sectionParsers_.end(); ++i) {
      i->isActive = false;
    }

    const int bufferSize = 65535;
    char buffer[bufferSize];
    int lineNo = 0;
    std::stack<std::string> sectionNameStack;
    // scan through the input stream and find section names
    while (input.getline(buffer, bufferSize)) {
      ++lineNo;

      std::string line = stripComments(Utils::trimLeftCopy(buffer));
      // a line begins with "//" is a comment line
      if (line.empty() ||
          (line.size() >= 2 && line[0] == '/' && line[1] == '/')) {
        continue;
      } else {
        StringTokenizer tokenizer(line);
        if (tokenizer.countTokens() < 2) {
          continue;
        } else {
          std::string keyword = tokenizer.nextToken();

          if (keyword == "begin") {
            std::string section = tokenizer.nextToken();
            sectionNameStack.push(section);

            i = std::find_if(sectionParsers_.begin(), sectionParsers_.end(),
                             SameSectionParserFunctor(section));
            if (i == sectionParsers_.end()) {
              snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                       "SectionParserManager Error: Can not find corresponding "
                       "section parser for %s\n",
                       section.c_str());
              painCave.isFatal = 1;
              simError();
            } else {
              if (i->isActive) {
                snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                         "SectionParserManager Error: Found multiple %s "
                         "sections\n",
                         section.c_str());
                painCave.isFatal = 1;
                simError();
              } else {
                i->isActive = true;
                i->lineNo   = lineNo;
                i->offset   = input.tellg();
              }
            }
          } else if (keyword == "end") {
            std::string section = tokenizer.nextToken();
            if (sectionNameStack.top() == section) {
              sectionNameStack.pop();
            } else {
              snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                       "SectionParserManager Error: begin %s "
                       "and end %s do not match at line %d\n",
                       sectionNameStack.top().c_str(), section.c_str(), lineNo);
              painCave.isFatal = 1;
              simError();
            }
          } else {
            continue;
          }
        }
      }
    }

    if (!sectionNameStack.empty()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "SectionParserManager Error: Stack is not empty.\n"
	       "\tCheck for matching begin / end lines.\n");
      painCave.isFatal = 1;
      simError();
    }

    // invoke parser
    for (i = sectionParsers_.begin(); i != sectionParsers_.end(); ++i) {
      if (i->isActive) {
        // C++ standard does not guarantee seekg resets EOF, in that
        // case, seekg will fail. It is always a good idea to call
        // clear() before seek
        input.clear();
        input.seekg(i->offset);
        (i->sectionParser)->parse(input, ff, i->lineNo);
        (i->sectionParser)->validateSection(ff);
      }
    }
  }

  void SectionParserManager::push_front(SectionParser* sp) {
    SectionParserManager::iterator i;
    i = findSectionParser(sp->getSectionName());
    if (i != sectionParsers_.end()) {
      std::cerr << sp->getSectionName() << " section parser already exists"
                << std::endl;
      return;
    }

    SectionParserContext context;

    if (sectionParsers_.empty()) {
      context.priority = beginPriority_;
    } else {
      context.priority = sectionParsers_.front().priority - priorityDifference_;
    }

    context.sectionParser = sp;
    context.lineNo        = 0;
    context.offset        = 0;
    context.isActive      = false;

    sectionParsers_.push_front(context);
  }

  void SectionParserManager::push_back(SectionParser* sp) {
    SectionParserManager::iterator i;
    i = findSectionParser(sp->getSectionName());
    if (i != sectionParsers_.end()) {
      std::cerr << sp->getSectionName() << " section parser already exists"
                << std::endl;
      return;
    }

    SectionParserContext context;
    if (sectionParsers_.empty()) {
      context.priority = beginPriority_;
    } else {
      context.priority = sectionParsers_.back().priority + priorityDifference_;
    }

    context.sectionParser = sp;
    context.lineNo        = 0;
    context.offset        = 0;
    context.isActive      = false;

    sectionParsers_.push_back(context);
  }

  void SectionParserManager::insert(SectionParser* sp, int priority) {
    SectionParserManager::iterator i;
    i = findSectionParser(sp->getSectionName());
    if (i != sectionParsers_.end()) {
      std::cerr << sp->getSectionName() << " section parser already exists"
                << std::endl;
    }

    SectionParserContext context;
    context.priority      = priority;
    context.sectionParser = sp;
    context.lineNo        = 0;
    context.offset        = 0;
    context.isActive      = false;

    if (sectionParsers_.empty()) {
      sectionParsers_.push_back(context);
    } else {
      for (i = sectionParsers_.begin(); i != sectionParsers_.end(); ++i) {
        if (i->priority == priority) {
          std::cerr << "Priority " << priority << " already used" << std::endl;
          return;
        } else if (i->priority > priority) {
          sectionParsers_.insert(i, context);
          break;
        }
      }
    }
  }

  SectionParserManager::iterator SectionParserManager::findSectionParser(
      const std::string& sectionName) {
    SectionParserManager::iterator i;
    for (i = sectionParsers_.begin(); i != sectionParsers_.end(); ++i) {
      if (i->sectionParser->getSectionName() == sectionName) { break; }
    }
    return i;
  }
}  // namespace OpenMD
