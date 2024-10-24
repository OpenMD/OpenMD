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

#ifndef OPENMD_UTILS_ACCUMULATOR_HPP
#define OPENMD_UTILS_ACCUMULATOR_HPP

#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include "math/SquareMatrix.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector.hpp"
#include "math/Vector3.hpp"
#include "nonbonded/NonBondedInteraction.hpp"
#include "utils/simError.h"

namespace OpenMD::Utils {

  template<typename T>
  class Accumulator {
    static_assert(std::is_default_constructible_v<T>,
                  "Accumulator type parameters must be default constructible.");

  public:
    void add(const T& val) {
      Count_++;

      Val_ = val;
      Total_ += val;
      Avg_ += (val - Avg_) / static_cast<RealType>(Count_);
      Avg2_ += (val * val - Avg2_) / static_cast<RealType>(Count_);

      if (Count_ <= 1) {
        Max_ = val;
        Min_ = val;
      } else {
        Max_ = val > Max_ ? val : Max_;
        Min_ = val < Min_ ? val : Min_;
      }
    }

    std::size_t getCount() const { return Count_; }

    T getLastValue() const { return Val_; }

    T getTotal() const {
      assert(Count_ != 0);
      return Total_;
    }

    RealType getMax() const {
      static_assert(std::is_arithmetic<T>::value,
                    "getMax() requires a RealType Accumulator.");
      assert(Count_ != 0);
      return Max_;
    }

    RealType getMin() const {
      static_assert(std::is_arithmetic<T>::value,
                    "getMin() requires a RealType Accumulator.");
      assert(Count_ != 0);
      return Min_;
    }

    RealType getAverage() const {
      assert(Count_ != 0);
      return Avg_;
    }

    RealType getVariance() const {
      assert(Count_ != 0);
      T var = (Avg2_ - Avg_ * Avg_);
      if (var < 0) var = 0;

      return var;
    }

    RealType getStdDev() const {
      assert(Count_ != 0);
      T sd = std::sqrt(this->getVariance());

      return sd;
    }

    RealType get95percentConfidenceInterval() const {
      assert(Count_ != 0);
      T ci =
          1.960 * this->getStdDev() / std::sqrt(static_cast<RealType>(Count_));

      return ci;
    }

  private:
    std::size_t Count_ {};
    T Val_ {}, Total_ {};
    RealType Max_ {}, Min_ {}, Avg_ {}, Avg2_ {};
  };

  // Specializations for commonly used Accumulator types
  template<>
  class Accumulator<std::vector<RealType>> {
  public:
    /* A flag specifying that a given bin is empty and should be ignored
      during calls to add() */
    constexpr static RealType BinEmptyFlag =
        std::numeric_limits<RealType>::max();

    void add(const std::vector<RealType>& val) {
      if (val.empty() || (val.size() != Avg_.size() && !Avg_.empty())) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Size of vector passed to add() did not "
                 "match the size of the StaticAccumulator.");
        painCave.isFatal = 1;
        simError();
      }

      Count_++;

      if (Avg_.empty()) {
        Val_.resize(val.size());
        Total_.resize(val.size());
        Avg_.resize(val.size());
        Avg2_.resize(val.size());
      }

      for (std::size_t i = 0; i < val.size(); i++) {
        /* If our placeholder, BinEmptyFlag, is passed to add(), we should
            not record data at the current index */
        if (val[i] == BinEmptyFlag) continue;
        Val_[i] = val[i];
        Total_[i] += val[i];
        Avg_[i] += (val[i] - Avg_[i]) / static_cast<RealType>(Count_);
        Avg2_[i] +=
            (val[i] * val[i] - Avg2_[i]) / static_cast<RealType>(Count_);
      }
    }

    std::size_t getCount() const { return Count_; }

    std::vector<RealType> getLastValue() const { return Val_; }

    std::vector<RealType> getTotal() const { return Total_; }

    std::vector<RealType> getAverage() const { return Avg_; }

    std::vector<RealType> getVariance() const {
      std::vector<RealType> var(Avg_.size());

      for (std::size_t i = 0; i < Avg_.size(); i++) {
        var[i] = (Avg2_[i] - Avg_[i] * Avg_[i]);
        if (var[i] < 0) var[i] = 0;
      }

      return var;
    }

    std::vector<RealType> getStdDev() const {
      std::vector<RealType> sd(Avg_.size());
      std::vector<RealType> variance = this->getVariance();

      for (std::size_t i = 0; i < variance.size(); i++) {
        sd[i] = std::sqrt(variance[i]);
      }

      return sd;
    }

    std::vector<RealType> get95percentConfidenceInterval() const {
      std::vector<RealType> ci(Avg_.size());
      std::vector<RealType> stdDev = this->getStdDev();

      for (std::size_t i = 0; i < stdDev.size(); i++) {
        ci[i] = 1.960 * stdDev[i] / std::sqrt(static_cast<RealType>(Count_));
      }

      return ci;
    }

  private:
    std::size_t Count_ {};
    std::vector<RealType> Val_ {}, Total_ {}, Avg_ {}, Avg2_ {};
  };

  template<unsigned int Dim>
  class Accumulator<Vector<RealType, Dim>> {
  public:
    void add(const Vector<RealType, Dim>& val) {
      Count_++;

      for (std::size_t i = 0; i < Dim; i++) {
        Val_[i] = val[i];
        Total_[i] += val[i];
        Avg_[i] += (val[i] - Avg_[i]) / static_cast<RealType>(Count_);
        Avg2_[i] +=
            (val[i] * val[i] - Avg2_[i]) / static_cast<RealType>(Count_);
      }
    }

    std::size_t getCount() const { return Count_; }

    Vector<RealType, Dim> getLastValue() const { return Val_; }

    Vector<RealType, Dim> getTotal() const {
      assert(Count_ != 0);

      return Total_;
    }

    Vector<RealType, Dim> getAverage() const {
      assert(Count_ != 0);

      return Avg_;
    }

    Vector<RealType, Dim> getVariance() const {
      assert(Count_ != 0);

      Vector<RealType, Dim> var {};

      for (std::size_t i = 0; i < Dim; i++) {
        var[i] = (Avg2_[i] - Avg_[i] * Avg_[i]);
        if (var[i] < 0) var[i] = 0;
      }

      return var;
    }

    Vector<RealType, Dim> getStdDev() const {
      assert(Count_ != 0);

      Vector<RealType, Dim> sd {};
      Vector<RealType, Dim> variance = this->getVariance();

      for (std::size_t i = 0; i < Dim; i++) {
        sd[i] = std::sqrt(variance[i]);
      }

      return sd;
    }

    Vector<RealType, Dim> get95percentConfidenceInterval() const {
      assert(Count_ != 0);

      Vector<RealType, Dim> ci {};
      Vector<RealType, Dim> stdDev = this->getStdDev();

      for (std::size_t i = 0; i < Dim; i++) {
        ci[i] = 1.960 * stdDev[i] / std::sqrt(static_cast<RealType>(Count_));
      }

      return ci;
    }

  private:
    std::size_t Count_ {};
    Vector<RealType, Dim> Val_ {}, Total_ {}, Avg_ {}, Avg2_ {};
  };

  template<>
  class Accumulator<Mat3x3d> {
  public:
    void add(const Mat3x3d& val) {
      Count_++;

      for (std::size_t i = 0; i < 3; i++) {
        for (std::size_t j = 0; j < 3; j++) {
          Val_(i, j) = val(i, j);
          Total_(i, j) += val(i, j);
          Avg_(i, j) +=
              (val(i, j) - Avg_(i, j)) / static_cast<RealType>(Count_);
          Avg2_(i, j) += (val(i, j) * val(i, j) - Avg2_(i, j)) /
                         static_cast<RealType>(Count_);
        }
      }
    }

    std::size_t getCount() const { return Count_; }

    Mat3x3d getLastValue() const { return Val_; }

    Mat3x3d getTotal() const {
      assert(Count_ != 0);

      return Total_;
    }

    Mat3x3d getAverage() const {
      assert(Count_ != 0);

      return Avg_;
    }

    Mat3x3d getVariance() const {
      assert(Count_ != 0);

      Mat3x3d var {};

      for (std::size_t i = 0; i < 3; i++) {
        for (std::size_t j = 0; j < 3; j++) {
          var(i, j) = (Avg2_(i, j) - Avg_(i, j) * Avg_(i, j));
          if (var(i, j) < 0) var(i, j) = 0;
        }
      }

      return var;
    }

    Mat3x3d getStdDev() const {
      assert(Count_ != 0);

      Mat3x3d sd {};
      Mat3x3d variance = this->getVariance();

      for (std::size_t i = 0; i < 3; i++) {
        for (std::size_t j = 0; j < 3; j++) {
          sd(i, j) = std::sqrt(variance(i, j));
        }
      }

      return sd;
    }

    Mat3x3d get95percentConfidenceInterval() const {
      assert(Count_ != 0);

      Mat3x3d ci {};
      Mat3x3d stdDev = this->getStdDev();

      for (std::size_t i = 0; i < 3; i++) {
        for (std::size_t j = 0; j < 3; j++) {
          ci(i, j) =
              1.960 * stdDev(i, j) / std::sqrt(static_cast<RealType>(Count_));
        }
      }

      return ci;
    }

  private:
    std::size_t Count_ {};
    Mat3x3d Val_ {}, Total_ {}, Avg_ {}, Avg2_ {};
  };

  // Type aliases for the most commonly used Accumulators
  using RealAccumulator      = Accumulator<RealType>;
  using StdVectorAccumulator = Accumulator<std::vector<RealType>>;
  using Vector3dAccumulator  = Accumulator<Vector<RealType, 3>>;
  using PotVecAccumulator =
      Accumulator<Vector<RealType, N_INTERACTION_FAMILIES>>;
  using Mat3x3dAccumulator = Accumulator<Mat3x3d>;
}  // namespace OpenMD::Utils

#endif  // OPENMD_UTILS_STATICACCUMULATOR_HPP
