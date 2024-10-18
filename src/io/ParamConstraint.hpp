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

#ifndef IO_PARAMCONSTRAINT_HPP
#define IO_PARAMCONSTRAINT_HPP

#include <sstream>

#include "utils/CaseConversion.hpp"
#include "utils/StringTokenizer.hpp"

namespace OpenMD {
  /**
   * This class allows to recognize constraint predicates, so that they can be
   * combined using composition operators. Every constraint predicate must be
   * derived from this class
   */
  template<typename Derived>
  struct ParamConstraintFacade {
    std::string getConstraintDescription() { return description_; }

  protected:
    std::string description_;
  };

  struct NotEmptyConstraint : public ParamConstraintFacade<NotEmptyConstraint> {
    NotEmptyConstraint() { description_ = "nonempty"; }
    bool operator()(const std::string& data) const { return !data.empty(); }
  };

  struct ZeroConstraint : public ParamConstraintFacade<ZeroConstraint> {
    ZeroConstraint() { this->description_ = "zero"; }
    template<typename DataType>
    bool operator()(DataType data) const {
      return data == 0;
    }
  };

  struct NonZeroConstraint : public ParamConstraintFacade<NonZeroConstraint> {
    NonZeroConstraint() { this->description_ = "nonzero"; }

    template<typename DataType>
    bool operator()(DataType data) const {
      return data != 0;
    }
  };

  struct PositiveConstraint : public ParamConstraintFacade<PositiveConstraint> {
    PositiveConstraint() { this->description_ = "positive"; }
    template<typename DataType>
    bool operator()(DataType data) const {
      return data > 0;
    }
  };

  struct NonPositiveConstraint :
      public ParamConstraintFacade<NonPositiveConstraint> {
    NonPositiveConstraint() { this->description_ = "nonpositive"; }
    template<typename DataType>
    bool operator()(DataType data) const {
      return data <= 0;
    }
  };

  struct NegativeConstraint : public ParamConstraintFacade<NegativeConstraint> {
    NegativeConstraint() { this->description_ = "negative"; }
    template<typename DataType>
    bool operator()(DataType data) const {
      return data < 0;
    }
  };

  struct NonNegativeConstraint :
      public ParamConstraintFacade<NonNegativeConstraint> {
    NonNegativeConstraint() { this->description_ = "nonnegative"; }
    template<typename DataType>
    bool operator()(DataType data) const {
      return data >= 0;
    }
  };

  struct EvenConstraint : public ParamConstraintFacade<EvenConstraint> {
    EvenConstraint() { this->description_ = "even"; }
    template<typename DataType>
    bool operator()(DataType data) const {
      return data % 2 == 0;
    }
  };

  template<typename T>
  struct LessThanConstraint :
      public ParamConstraintFacade<LessThanConstraint<T>> {
    LessThanConstraint(T rhs) : rhs_(rhs) {
      std::stringstream iss;
      iss << "less than " << rhs;
      this->description_ = iss.str();
    }
    template<typename DataType>
    bool operator()(DataType data) const {
      return data < rhs_;
    }

  private:
    T rhs_;
  };

  template<typename T>
  struct LessThanOrEqualToConstraint :
      public ParamConstraintFacade<LessThanOrEqualToConstraint<T>> {
    LessThanOrEqualToConstraint(T rhs) : rhs_(rhs) {
      std::stringstream iss;
      iss << "less than or equal to" << rhs;
      this->description_ = iss.str();
    }

    template<typename DataType>
    bool operator()(DataType data) const {
      return data <= rhs_;
    }

  private:
    T rhs_;
  };

  template<typename T>
  struct EqualConstraint : public ParamConstraintFacade<EqualConstraint<T>> {
    EqualConstraint(T rhs) : rhs_(rhs) {
      std::stringstream iss;
      iss << "equal to" << rhs;
      this->description_ = iss.str();
    }
    template<typename DataType>
    bool operator()(DataType data) const {
      return data == rhs_;
    }

  private:
    T rhs_;
  };

  struct EqualIgnoreCaseConstraint :
      public ParamConstraintFacade<EqualIgnoreCaseConstraint> {
    EqualIgnoreCaseConstraint(std::string rhs) :
        rhs_(OpenMD::toUpperCopy(rhs)) {
      std::stringstream iss;
      iss << "equal to (case insensitive) " << rhs;
      this->description_ = iss.str();
    }

    bool operator()(std::string data) const {
      return OpenMD::toUpperCopy(data) == rhs_;
    }

  private:
    std::string rhs_;
  };

  struct ContainsConstraint :
      public ParamConstraintFacade<EqualIgnoreCaseConstraint> {
    ContainsConstraint(std::string rhs) : rhs_(OpenMD::toUpperCopy(rhs)) {
      std::stringstream iss;
      iss << "contains " << rhs;
      this->description_ = iss.str();
    }

    bool operator()(std::string data) const {
      OpenMD::StringTokenizer tokenizer(OpenMD::toUpperCopy(data),
                                        " ,;|\t\n\r");
      while (tokenizer.hasMoreTokens()) {
        if (tokenizer.nextToken() == rhs_) { return true; }
      }

      return false;
    }

  private:
    std::string rhs_;
  };

  template<typename T>
  struct GreaterThanConstraint :
      public ParamConstraintFacade<GreaterThanConstraint<T>> {
    GreaterThanConstraint(T rhs) : rhs_(rhs) {
      std::stringstream iss;
      iss << "greater than" << rhs;
      this->description_ = iss.str();
    }
    template<typename DataType>
    bool operator()(DataType data) const {
      return data > rhs_;
    }

  private:
    T rhs_;
  };

  template<typename T>
  struct GreaterThanOrEqualTo :
      public ParamConstraintFacade<GreaterThanOrEqualTo<T>> {
    GreaterThanOrEqualTo(T rhs) : rhs_(rhs) {
      std::stringstream iss;
      iss << "greater than or equal to" << rhs;
      this->description_ = iss.str();
    }
    template<typename DataType>
    bool operator()(DataType data) const {
      return data >= rhs_;
    }

  private:
    T rhs_;
  };

  // class_and composition predicate
  template<typename Cons1T, typename Cons2T>
  struct AndParamConstraint :
      public ParamConstraintFacade<AndParamConstraint<Cons1T, Cons2T>> {
  public:
    AndParamConstraint(Cons1T cons1, Cons2T cons2) :
        cons1_(cons1), cons2_(cons2) {
      std::stringstream iss;
      iss << "(" << cons1_.getConstraintDescription() << " and "
          << cons2_.getConstraintDescription() << ")";
      this->description_ = iss.str();
    }

    template<typename DataType>
    bool operator()(DataType data) const {
      return cons1_(data) && cons2_(data);
    }

  private:
    Cons1T cons1_;
    Cons2T cons2_;
  };

  template<typename Cons1T, typename Cons2T>
  struct OrParamConstraint :
      public ParamConstraintFacade<OrParamConstraint<Cons1T, Cons2T>> {
  public:
    OrParamConstraint(Cons1T cons1, Cons2T cons2) :
        cons1_(cons1), cons2_(cons2) {
      std::stringstream iss;
      iss << cons1_.getConstraintDescription() << " or "
          << cons2_.getConstraintDescription() << "";
      this->description_ = iss.str();
    }

    template<typename DataType>
    bool operator()(DataType data) const {
      return cons1_(data) || cons2_(data);
    }

  private:
    Cons1T cons1_;
    Cons2T cons2_;
  };

  template<typename ConsT>
  struct NotParamConstraint :
      public ParamConstraintFacade<NotParamConstraint<ConsT>> {
  public:
    NotParamConstraint(ConsT cons) : cons_(cons) {
      std::stringstream iss;
      iss << "(not" << cons_.getConstraintDescription() << ")";
      this->description_ = iss.str();
    }

    template<typename DataType>
    bool operator()(DataType data) const {
      return !cons_(data);
    }

  private:
    ConsT cons_;
  };

  template<typename Cons1T, typename Cons2T>
  inline AndParamConstraint<Cons1T, Cons2T> operator&&(
      const ParamConstraintFacade<Cons1T>& cons1,
      const ParamConstraintFacade<Cons2T>& cons2) {
    return AndParamConstraint<Cons1T, Cons2T>(
        *static_cast<const Cons1T*>(&cons1),
        *static_cast<const Cons2T*>(&cons2));
  }

  template<typename Cons1T, typename Cons2T>
  inline OrParamConstraint<Cons1T, Cons2T> operator||(
      const ParamConstraintFacade<Cons1T>& cons1,
      const ParamConstraintFacade<Cons2T>& cons2) {
    return OrParamConstraint<Cons1T, Cons2T>(
        *static_cast<const Cons1T*>(&cons1),
        *static_cast<const Cons2T*>(&cons2));
  }

  template<typename ConsT>
  inline NotParamConstraint<ConsT> operator!(
      const ParamConstraintFacade<ConsT>& cons) {
    return NotParamConstraint<ConsT>(*static_cast<const ConsT*>(&cons));
  }

  NotEmptyConstraint isNotEmpty();
  ZeroConstraint isZero();

  ParamConstraintFacade<NonZeroConstraint> isNonZero();
  PositiveConstraint isPositive();
  NonPositiveConstraint isNonPositive();

  NegativeConstraint isNegative();

  NonNegativeConstraint isNonNegative();
  EvenConstraint isEven();

  template<typename T>
  inline LessThanConstraint<T> isLessThan(T& v) {
    return LessThanConstraint<T>(v);
  }

  template<typename T>
  inline LessThanOrEqualToConstraint<T> isLessThanOrEqualTo(T& v) {
    return LessThanOrEqualToConstraint<T>(v);
  }

  template<typename T>
  inline EqualConstraint<T> isEqual(T& v) {
    return EqualConstraint<T>(v);
  }

  template<typename T>
  inline GreaterThanConstraint<T> isGreaterThan(T& v) {
    return GreaterThanConstraint<T>(v);
  }

  template<typename T>
  inline GreaterThanOrEqualTo<T> isGreaterThanOrEqualTo(T& v) {
    return GreaterThanOrEqualTo<T>(v);
  }

  EqualIgnoreCaseConstraint isEqualIgnoreCase(std::string str);
}  // namespace OpenMD

#endif
