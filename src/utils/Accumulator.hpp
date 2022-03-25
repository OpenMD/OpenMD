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

#ifndef UTILS_ACCUMULATOR_HPP
#define UTILS_ACCUMULATOR_HPP

#include <cassert>
#include <cmath>

#include "math/Vector3.hpp"
#include "nonbonded/NonBondedInteraction.hpp"

namespace OpenMD {

  class BaseAccumulator {
  public:
    virtual void clear() = 0;
    /**
     * get the number of accumulated values
     */
    virtual size_t count() { return Count_; }
    virtual ~BaseAccumulator() {};

  protected:
    size_t Count_;
  };

  /**
   * Basic Accumulator class for numbers.
   */
  class Accumulator : public BaseAccumulator {
    using ElementType = RealType;
    using ResultType  = RealType;

  public:
    Accumulator() : BaseAccumulator() { this->clear(); }

    ~Accumulator() {};

    /**
     * Accumulate another value
     */
    virtual void add(ElementType const& val) {
      Count_++;
      Avg_ += (val - Avg_) / Count_;
      Avg2_ += (val * val - Avg2_) / Count_;
      Val_ = val;
      Total_ += val;
      if (Count_ <= 1) {
        Max_ = val;
        Min_ = val;
      } else {
        Max_ = val > Max_ ? val : Max_;
        Min_ = val < Min_ ? val : Min_;
      }
    }

    /**
     * reset the Accumulator to the empty state
     */
    void clear() {
      Count_ = 0;
      Avg_   = 0;
      Avg2_  = 0;
      Total_ = 0;
      Val_   = 0;
    }

    /**
     * return the most recently added value
     */
    void getLastValue(ElementType& ret) {
      ret = Val_;
      return;
    }

    /**
     * compute the Mean
     */
    void getAverage(ResultType& ret) {
      assert(Count_ != 0);
      ret = Avg_;
      return;
    }

    /**
     * return the Total accumulated sum
     */
    void getTotal(ResultType& ret) {
      assert(Count_ != 0);
      ret = Total_;
      return;
    }

    /**
     * compute the Variance
     */
    void getVariance(ResultType& ret) {
      assert(Count_ != 0);
      ret = (Avg2_ - Avg_ * Avg_);
      if (ret < 0) ret = 0;
      return;
    }

    /**
     * compute error of average value
     */
    void getStdDev(ResultType& ret) {
      assert(Count_ != 0);
      RealType var;
      this->getVariance(var);
      ret = sqrt(var);
      return;
    }

    /**
     * return the largest value
     */
    void getMax(ElementType& ret) {
      assert(Count_ != 0);
      ret = Max_;
      return;
    }

    /**
     * return the smallest value
     */
    void getMin(ElementType& ret) {
      assert(Count_ != 0);
      ret = Max_;
      return;
    }

    /**
     * return the 95% confidence interval:
     *
     * That is returns c, such that we have 95% confidence that the
     * true mean is within 2c of the Average (x):
     *
     *   x - c <= true mean <= x + c
     *
     */
    void get95percentConfidenceInterval(ResultType& ret) {
      assert(Count_ != 0);
      RealType sd;
      this->getStdDev(sd);
      ret = 1.960 * sd / sqrt(RealType(Count_));
      return;
    }

  private:
    ElementType Val_;
    ResultType Avg_;
    ResultType Avg2_;
    ResultType Total_;
    ElementType Min_;
    ElementType Max_;
  };

  class VectorAccumulator : public BaseAccumulator {
    using ElementType = Vector3d;
    using ResultType  = Vector3d;

  public:
    VectorAccumulator() : BaseAccumulator() { this->clear(); }

    /**
     * Accumulate another value
     */
    void add(ElementType const& val) {
      Count_++;
      RealType len(0.0);
      for (unsigned int i = 0; i < 3; i++) {
        Avg_[i] += (val[i] - Avg_[i]) / Count_;
        Avg2_[i] += (val[i] * val[i] - Avg2_[i]) / Count_;
        Val_[i] = val[i];
        Total_[i] += val[i];
        len += val[i] * val[i];
      }
      len = sqrt(len);
      AvgLen_ += (len - AvgLen_) / Count_;
      AvgLen2_ += (len * len - AvgLen2_) / Count_;

      if (Count_ <= 1) {
        Max_ = len;
        Min_ = len;
      } else {
        Max_ = len > Max_ ? len : Max_;
        Min_ = len < Min_ ? len : Min_;
      }
    }

    /**
     * reset the Accumulator to the empty state
     */
    void clear() {
      Count_   = 0;
      Avg_     = V3Zero;
      Avg2_    = V3Zero;
      Total_   = V3Zero;
      Val_     = V3Zero;
      AvgLen_  = 0;
      AvgLen2_ = 0;
    }

    /**
     * return the most recently added value
     */
    void getLastValue(ElementType& ret) {
      ret = Val_;
      return;
    }

    /**
     * compute the Mean
     */
    void getAverage(ResultType& ret) {
      assert(Count_ != 0);
      ret = Avg_;
      return;
    }

    /**
     * return the Total accumulated sum
     */
    void getTotal(ResultType& ret) {
      assert(Count_ != 0);
      ret = Total_;
      return;
    }

    /**
     * compute the Variance
     */
    void getVariance(ResultType& ret) {
      assert(Count_ != 0);
      for (unsigned int i = 0; i < 3; i++) {
        ret[i] = (Avg2_[i] - Avg_[i] * Avg_[i]);
        if (ret[i] < 0) ret[i] = 0;
      }
      return;
    }

    /**
     * compute error of average value
     */
    void getStdDev(ResultType& ret) {
      assert(Count_ != 0);
      ResultType var;
      this->getVariance(var);
      ret[0] = sqrt(var[0]);
      ret[1] = sqrt(var[1]);
      ret[2] = sqrt(var[2]);
      return;
    }

    /**
     * return the 95% confidence interval:
     *
     * That is returns c, such that we have 95% confidence that the
     * true mean is within 2c of the Average (x):
     *
     *   x - c <= true mean <= x + c
     *
     */
    void get95percentConfidenceInterval(ResultType& ret) {
      assert(Count_ != 0);
      ResultType sd;
      this->getStdDev(sd);
      ret[0] = 1.960 * sd[0] / sqrt(RealType(Count_));
      ret[1] = 1.960 * sd[1] / sqrt(RealType(Count_));
      ret[2] = 1.960 * sd[2] / sqrt(RealType(Count_));
      return;
    }

    /**
     * return the largest length
     */
    void getMaxLength(RealType& ret) {
      assert(Count_ != 0);
      ret = Max_;
      return;
    }

    /**
     * return the smallest length
     */
    void getMinLength(RealType& ret) {
      assert(Count_ != 0);
      ret = Min_;
      return;
    }

    /**
     * return the largest length
     */
    void getAverageLength(RealType& ret) {
      assert(Count_ != 0);
      ret = AvgLen_;
      return;
    }

    /**
     * compute the Variance of the length
     */
    void getLengthVariance(RealType& ret) {
      assert(Count_ != 0);
      ret = (AvgLen2_ - AvgLen_ * AvgLen_);
      if (ret < 0) ret = 0;
      return;
    }

    /**
     * compute error of average value
     */
    void getLengthStdDev(RealType& ret) {
      assert(Count_ != 0);
      RealType var;
      this->getLengthVariance(var);
      ret = sqrt(var);
      return;
    }

    /**
     * return the 95% confidence interval:
     *
     * That is returns c, such that we have 95% confidence that the
     * true mean is within 2c of the Average (x):
     *
     *   x - c <= true mean <= x + c
     *
     */
    void getLength95percentConfidenceInterval(ResultType& ret) {
      assert(Count_ != 0);
      RealType sd;
      this->getLengthStdDev(sd);
      ret = 1.960 * sd / sqrt(RealType(Count_));
      return;
    }

  private:
    ResultType Val_;
    ResultType Avg_;
    ResultType Avg2_;
    ResultType Total_;
    RealType AvgLen_;
    RealType AvgLen2_;
    RealType Min_;
    RealType Max_;
  };

  class PotVecAccumulator : public BaseAccumulator {
    using ElementType = potVec;
    using ResultType  = potVec;

  public:
    PotVecAccumulator() : BaseAccumulator() { this->clear(); }

    /**
     * Accumulate another value
     */
    void add(ElementType const& val) {
      Count_++;
      RealType len(0.0);
      for (unsigned int i = 0; i < N_INTERACTION_FAMILIES; i++) {
        Avg_[i] += (val[i] - Avg_[i]) / Count_;
        Avg2_[i] += (val[i] * val[i] - Avg2_[i]) / Count_;
        Val_[i] = val[i];
        Total_[i] += val[i];
        len += val[i] * val[i];
      }
      len = sqrt(len);
      AvgLen_ += (len - AvgLen_) / Count_;
      AvgLen2_ += (len * len - AvgLen2_) / Count_;

      if (Count_ <= 1) {
        Max_ = len;
        Min_ = len;
      } else {
        Max_ = len > Max_ ? len : Max_;
        Min_ = len < Min_ ? len : Min_;
      }
    }

    /**
     * reset the Accumulator to the empty state
     */
    void clear() {
      Count_ = 0;
      const Vector<RealType, N_INTERACTION_FAMILIES> potVecZero(0.0);
      Avg_     = potVecZero;
      Avg2_    = potVecZero;
      Val_     = potVecZero;
      Total_   = potVecZero;
      AvgLen_  = 0;
      AvgLen2_ = 0;
    }

    /**
     * return the most recently added value
     */
    void getLastValue(ElementType& ret) {
      ret = Val_;
      return;
    }

    /**
     * compute the Mean
     */
    void getAverage(ResultType& ret) {
      assert(Count_ != 0);
      ret = Avg_;
      return;
    }

    /**
     * return the Total accumulated sum
     */
    void getTotal(ResultType& ret) {
      assert(Count_ != 0);
      ret = Total_;
      return;
    }
    /**
     * compute the Variance
     */
    void getVariance(ResultType& ret) {
      assert(Count_ != 0);
      for (unsigned int i = 0; i < N_INTERACTION_FAMILIES; i++) {
        ret[i] = (Avg2_[i] - Avg_[i] * Avg_[i]);
        if (ret[i] < 0) ret[i] = 0;
      }
      return;
    }

    /**
     * compute error of average value
     */
    void getStdDev(ResultType& ret) {
      assert(Count_ != 0);
      ResultType var;
      this->getVariance(var);
      for (unsigned int i = 0; i < N_INTERACTION_FAMILIES; i++) {
        ret[i] = sqrt(var[i]);
      }
      return;
    }

    /**
     * return the 95% confidence interval:
     *
     * That is returns c, such that we have 95% confidence that the
     * true mean is within 2c of the Average (x):
     *
     *   x - c <= true mean <= x + c
     *
     */
    void get95percentConfidenceInterval(ResultType& ret) {
      assert(Count_ != 0);
      ResultType sd;
      this->getStdDev(sd);
      for (unsigned int i = 0; i < N_INTERACTION_FAMILIES; i++) {
        ret[i] = 1.960 * sd[i] / sqrt(RealType(Count_));
      }
      return;
    }

    /**
     * return the largest length
     */
    void getMaxLength(RealType& ret) {
      assert(Count_ != 0);
      ret = Max_;
      return;
    }

    /**
     * return the smallest length
     */
    void getMinLength(RealType& ret) {
      assert(Count_ != 0);
      ret = Min_;
      return;
    }

    /**
     * return the largest length
     */
    void getAverageLength(RealType& ret) {
      assert(Count_ != 0);
      ret = AvgLen_;
      return;
    }

    /**
     * compute the Variance of the length
     */
    void getLengthVariance(RealType& ret) {
      assert(Count_ != 0);
      ret = (AvgLen2_ - AvgLen_ * AvgLen_);
      if (ret < 0) ret = 0;
      return;
    }

    /**
     * compute error of average value
     */
    void getLengthStdDev(RealType& ret) {
      assert(Count_ != 0);
      RealType var;
      this->getLengthVariance(var);
      ret = sqrt(var);
      return;
    }

    /**
     * return the 95% confidence interval:
     *
     * That is returns c, such that we have 95% confidence that the
     * true mean is within 2c of the Average (x):
     *
     *   x - c <= true mean <= x + c
     *
     */
    void getLength95percentConfidenceInterval(ResultType& ret) {
      assert(Count_ != 0);
      RealType sd;
      this->getLengthStdDev(sd);
      ret = 1.960 * sd / sqrt(RealType(Count_));
      return;
    }

  private:
    ResultType Val_;
    ResultType Avg_;
    ResultType Avg2_;
    ResultType Total_;
    RealType AvgLen_;
    RealType AvgLen2_;
    RealType Min_;
    RealType Max_;
  };

  class MatrixAccumulator : public BaseAccumulator {
    using ElementType = Mat3x3d;
    using ResultType  = Mat3x3d;

  public:
    MatrixAccumulator() : BaseAccumulator() { this->clear(); }

    /**
     * Accumulate another value
     */
    void add(ElementType const& val) {
      Count_++;
      for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
          Avg_(i, j) += (val(i, j) - Avg_(i, j)) / Count_;
          Avg2_(i, j) += (val(i, j) * val(i, j) - Avg2_(i, j)) / Count_;
          Val_(i, j) = val(i, j);
          Total_(i, j) += val(i, j);
        }
      }
    }

    /**
     * reset the Accumulator to the empty state
     */
    void clear() {
      Count_ = 0;
      Avg_ *= 0.0;
      Avg2_ *= 0.0;
      Val_ *= 0.0;
      Total_ *= 0.0;
    }

    /**
     * return the most recently added value
     */
    void getLastValue(ElementType& ret) {
      ret = Val_;
      return;
    }

    /**
     * compute the Mean
     */
    void getAverage(ResultType& ret) {
      assert(Count_ != 0);
      ret = Avg_;
      return;
    }

    /**
     * return the Total accumulated sum
     */
    void getTotal(ResultType& ret) {
      assert(Count_ != 0);
      ret = Total_;
      return;
    }

    /**
     * compute the Variance
     */
    void getVariance(ResultType& ret) {
      assert(Count_ != 0);
      for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
          ret(i, j) = (Avg2_(i, j) - Avg_(i, j) * Avg_(i, j));
          if (ret(i, j) < 0) ret(i, j) = 0;
        }
      }
      return;
    }

    /**
     * compute error of average value
     */
    void getStdDev(ResultType& ret) {
      assert(Count_ != 0);
      Mat3x3d var;
      this->getVariance(var);
      for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
          ret(i, j) = sqrt(var(i, j));
        }
      }
      return;
    }

    /**
     * return the 95% confidence interval:
     *
     * That is returns c, such that we have 95% confidence that the
     * true mean is within 2c of the Average (x):
     *
     *   x - c <= true mean <= x + c
     *
     */
    void get95percentConfidenceInterval(ResultType& ret) {
      assert(Count_ != 0);
      Mat3x3d sd;
      this->getStdDev(sd);
      for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
          ret(i, j) = 1.960 * sd(i, j) / sqrt(RealType(Count_));
        }
      }
      return;
    }

  private:
    ElementType Val_;
    ResultType Avg_;
    ResultType Avg2_;
    ResultType Total_;
  };
}  // namespace OpenMD

#endif
