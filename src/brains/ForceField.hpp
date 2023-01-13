/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

/**
 * @file ForceField.hpp
 * @author tlin
 * @date 11/04/2004
 * @version 1.0
 */

#ifndef USETHEFORCE_FORCEFIELD_HPP
#define USETHEFORCE_FORCEFIELD_HPP

#include <config.h>

#include <string>
#include <utility>
#include <vector>

#include "io/ForceFieldOptions.hpp"
#include "io/SectionParserManager.hpp"
#include "io/ifstrstream.hpp"
#include "types/AtomType.hpp"
#include "types/BendType.hpp"
#include "types/BondType.hpp"
#include "types/InversionType.hpp"
#include "types/NonBondedInteractionType.hpp"
#include "types/TorsionType.hpp"
#include "utils/TypeContainer.hpp"

namespace OpenMD {

  class ForceField {
  public:
    using AtomTypeContainer      = TypeContainer<AtomType, 1>;
    using BondTypeContainer      = TypeContainer<BondType, 2>;
    using BendTypeContainer      = TypeContainer<BendType, 3>;
    using TorsionTypeContainer   = TypeContainer<TorsionType, 4>;
    using InversionTypeContainer = TypeContainer<InversionType, 4>;
    using NonBondedInteractionTypeContainer =
        TypeContainer<NonBondedInteractionType, 2>;

    ForceField(std::string ffName);

    virtual ~ForceField() = default;

    std::string getForceFieldFileName() { return forceFieldFileName_; }

    void setForceFieldFileName(const std::string& filename) {
      forceFieldFileName_ = filename;
    }

    virtual void parse(const std::string& filename);

    AtomType* getAtomType(const std::string& at);
    AtomType* getAtomType(int ident);
    BondType* getBondType(const std::string& at1, const std::string& at2);
    BendType* getBendType(const std::string& at1, const std::string& at2,
                          const std::string& at3);
    TorsionType* getTorsionType(const std::string& at1, const std::string& at2,
                                const std::string& at3, const std::string& at4);
    InversionType* getInversionType(const std::string& at1,
                                    const std::string& at2,
                                    const std::string& at3,
                                    const std::string& at4);
    NonBondedInteractionType* getNonBondedInteractionType(
        const std::string& at1, const std::string& at2);

    BondType* getExactBondType(const std::string& at1, const std::string& at2);
    BendType* getExactBendType(const std::string& at1, const std::string& at2,
                               const std::string& at3);
    TorsionType* getExactTorsionType(const std::string& at1,
                                     const std::string& at2,
                                     const std::string& at3,
                                     const std::string& at4);
    InversionType* getExactInversionType(const std::string& at1,
                                         const std::string& at2,
                                         const std::string& at3,
                                         const std::string& at4);
    NonBondedInteractionType* getExactNonBondedInteractionType(
        const std::string& at1, const std::string& at2);

    // avoid make virtual function public
    // Herb Sutter and Andrei Alexandrescu, C++ coding Standards,
    // Addision-Wesley
    virtual RealType getRcutFromAtomType(AtomType* at);

    std::string getWildCard() { return wildCardAtomTypeName_; }

    void setWildCard(const std::string& wildCard) {
      wildCardAtomTypeName_ = wildCard;
    }

    size_t getNAtomType() { return atomTypeCont_.size(); }

    AtomTypeContainer* getAtomTypes() { return &atomTypeCont_; }

    NonBondedInteractionTypeContainer* getNonBondedInteractionTypes() {
      return &nonBondedInteractionTypeCont_;
    }

    bool addAtomType(const std::string& at, AtomType* atomType);

    bool replaceAtomType(const std::string& at, AtomType* atomType);

    bool addBondType(const std::string& at1, const std::string& at2,
                     BondType* bondType);

    bool addBendType(const std::string& at1, const std::string& at2,
                     const std::string& at3, BendType* bendType);

    bool addTorsionType(const std::string& at1, const std::string& at2,
                        const std::string& at3, const std::string& at4,
                        TorsionType* torsionType);

    bool addInversionType(const std::string& at1, const std::string& at2,
                          const std::string& at3, const std::string& at4,
                          InversionType* inversionType);

    bool addNonBondedInteractionType(const std::string& at1,
                                     const std::string& at2,
                                     NonBondedInteractionType* nbiType);

    ifstrstream* openForceFieldFile(const std::string& filename);

    ForceFieldOptions& getForceFieldOptions() { return forceFieldOptions_; }

  protected:
    AtomTypeContainer atomTypeCont_;
    BondTypeContainer bondTypeCont_;
    BendTypeContainer bendTypeCont_;
    TorsionTypeContainer torsionTypeCont_;
    InversionTypeContainer inversionTypeCont_;
    NonBondedInteractionTypeContainer nonBondedInteractionTypeCont_;
    ForceFieldOptions forceFieldOptions_;
    std::map<int, std::string> atypeIdentToName;
    SectionParserManager spMan_;

  private:
    std::string ffPath_;
    std::string wildCardAtomTypeName_;
    std::string forceFieldFileName_;
  };
}  // namespace OpenMD

#endif
