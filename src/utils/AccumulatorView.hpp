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

#ifndef OPENMD_UTILS_ACCUMULATORVIEW_HPP
#define OPENMD_UTILS_ACCUMULATORVIEW_HPP

#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <typeinfo>
#include <vector>

#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "nonbonded/NonBondedInteraction.hpp"
#include "utils/Accumulator.hpp"
#include "utils/BaseAccumulator.hpp"
#include "utils/simError.h"

namespace OpenMD::Utils {

  template<typename T>
  class AccumulatorView;

  template<>
  class AccumulatorView<RealAccumulator> :
      public BaseAccumulator,
      private RealAccumulator {
  public:
    void add(RealType val) override { RealAccumulator::add(val); }

    // Other overrides for invalid entries
    void add(const std::vector<RealType>&) override {
      accumulatorFunctionCallMismatch();
    }
    void add(const Vector3d&) override { accumulatorFunctionCallMismatch(); }
    void add(const potVec&) override { accumulatorFunctionCallMismatch(); }
    void add(const Mat3x3d&) override { accumulatorFunctionCallMismatch(); }

    void writeData(std::ostream& stream, const std::string& errorMessage,
                   DataHandling dataHandling) const override {
      RealType dat {};
      std::size_t count = RealAccumulator::getCount();

      switch (dataHandling) {
      case DataHandling::Average:
        dat = RealAccumulator::getAverage();
        break;
      case DataHandling::Last:
        dat = RealAccumulator::getLastValue();
        break;
      case DataHandling::Max:
        dat = RealAccumulator::getMax();
        break;
      case DataHandling::Min:
        dat = RealAccumulator::getMin();
        break;
      case DataHandling::Total:
        dat = RealAccumulator::getTotal();
        break;
      default:
        break;
      }

      if (std::isinf(dat) || std::isnan(dat)) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "%s",
                 errorMessage.c_str());
        painCave.isFatal = 1;
        simError();
      } else {
        if (count == 0)
          stream << "\t";
        else
          stream << "\t" << dat;
      }
    }

    void writeErrorBars(std::ostream& stream, const std::string& errorMessage,
                        ErrorHandling errorHandling) const override {
      RealType err {};
      std::size_t count = RealAccumulator::getCount();

      switch (errorHandling) {
      case ErrorHandling::CI95:
        err = RealAccumulator::get95percentConfidenceInterval();
        break;
      case ErrorHandling::StdDev:
        err = RealAccumulator::getStdDev();
        break;
      case ErrorHandling::Variance:
        err = RealAccumulator::getVariance();
        break;
      default:
        break;
      }

      if (count == 0 && errorHandling == ErrorHandling::CI95) {
        stream << "\t";
      } else {
        if (std::isinf(err) || std::isnan(err)) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "%s",
                   errorMessage.c_str());
          painCave.isFatal = 1;
          simError();
        } else {
          if (count == 0)
            stream << "\t";
          else
            stream << "\t" << err;
        }
      }
    }

    std::type_index getType() const override { return typeid(RealType); }

    std::size_t getCount() const override {
      return RealAccumulator::getCount();
    }
  };

  template<>
  class AccumulatorView<StdVectorAccumulator> :
      public BaseAccumulator,
      private StdVectorAccumulator {
  public:
    void add(const std::vector<RealType>& val) override {
      StdVectorAccumulator::add(val);
    }

    // Other overrides for invalid entries
    void add(RealType val) override { accumulatorFunctionCallMismatch(); }
    void add(const Vector3d&) override { accumulatorFunctionCallMismatch(); }
    void add(const potVec&) override { accumulatorFunctionCallMismatch(); }
    void add(const Mat3x3d&) override { accumulatorFunctionCallMismatch(); }

    void writeData(std::ostream& stream, const std::string& errorMessage,
                   DataHandling dataHandling) const override {
      std::vector<RealType> dat;
      std::size_t count = StdVectorAccumulator::getCount();

      switch (dataHandling) {
      case DataHandling::Average:
        dat = StdVectorAccumulator::getAverage();
        break;
      case DataHandling::Last:
        dat = StdVectorAccumulator::getLastValue();
        break;
      case DataHandling::Total:
        dat = StdVectorAccumulator::getTotal();
        break;
      case DataHandling::Max:
      case DataHandling::Min:
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Max and Min functions are not defined for a "
                 "std::vector Accumulator.");
        painCave.isFatal = 1;
        simError();
        break;
      default:
        break;
      }

      for (int i = 0; i < StdVectorAccumulator::getAverage().size(); i++) {
        if (std::isinf(dat[i]) || std::isnan(dat[i])) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "%s",
                   errorMessage.c_str());
          painCave.isFatal = 1;
          simError();
        } else {
          if (count == 0)
            stream << '\t';
          else
            stream << '\t' << dat[i];
        }
      }
    }

    void writeErrorBars(std::ostream& stream, const std::string& errorMessage,
                        ErrorHandling errorHandling) const override {
      std::vector<RealType> err;
      std::size_t count = StdVectorAccumulator::getCount();

      switch (errorHandling) {
      case ErrorHandling::CI95:
        err = StdVectorAccumulator::get95percentConfidenceInterval();
        break;
      case ErrorHandling::StdDev:
        err = StdVectorAccumulator::getStdDev();
        break;
      case ErrorHandling::Variance:
        err = StdVectorAccumulator::getVariance();
        break;
      default:
        break;
      }

      for (int i = 0; i < StdVectorAccumulator::getAverage().size(); i++) {
        if (count == 0 && errorHandling == ErrorHandling::CI95) {
          stream << "\t";
        } else {
          if (std::isinf(err[i]) || std::isnan(err[i])) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "%s",
                     errorMessage.c_str());
            painCave.isFatal = 1;
            simError();
          } else {
            if (count == 0)
              stream << "\t";
            else
              stream << "\t" << err[i];
          }
        }
      }
    }

    std::type_index getType() const override {
      return typeid(std::vector<RealType>);
    }

    std::size_t getCount() const override {
      std::size_t count = StdVectorAccumulator::getCount();

      return count;
    }
  };

  template<>
  class AccumulatorView<Vector3dAccumulator> :
      public BaseAccumulator,
      private Vector3dAccumulator {
  public:
    void add(const Vector3d& val) override { Vector3dAccumulator::add(val); }

    // Other overrides for invalid entries
    void add(RealType val) override { accumulatorFunctionCallMismatch(); }
    void add(const std::vector<RealType>&) override {
      accumulatorFunctionCallMismatch();
    }
    void add(const potVec&) override { accumulatorFunctionCallMismatch(); }
    void add(const Mat3x3d&) override { accumulatorFunctionCallMismatch(); }

    void writeData(std::ostream& stream, const std::string& errorMessage,
                   DataHandling dataHandling) const override {
      Vector3d dat;
      std::size_t count = Vector3dAccumulator::getCount();

      switch (dataHandling) {
      case DataHandling::Average:
        dat = Vector3dAccumulator::getAverage();
        break;
      case DataHandling::Last:
        dat = Vector3dAccumulator::getLastValue();
        break;
      case DataHandling::Total:
        dat = Vector3dAccumulator::getTotal();
        break;
      case DataHandling::Max:
      case DataHandling::Min:
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "%s",
                 "Max and Min functions are not defined for a "
                 "std::vector Accumulator.");
        painCave.isFatal = 1;
        simError();
        break;
      default:
        break;
      }

      for (int i = 0; i < Vector3dAccumulator::getAverage().size(); i++) {
        if (std::isinf(dat[i]) || std::isnan(dat[i])) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "%s",
                   errorMessage.c_str());
          painCave.isFatal = 1;
          simError();
        } else {
          if (count == 0)
            stream << "\t";
          else
            stream << "\t" << dat[i];
        }
      }
    }

    void writeErrorBars(std::ostream& stream, const std::string& errorMessage,
                        ErrorHandling errorHandling) const override {
      Vector3d err;
      std::size_t count = Vector3dAccumulator::getCount();

      switch (errorHandling) {
      case ErrorHandling::CI95:
        err = Vector3dAccumulator::get95percentConfidenceInterval();
        break;
      case ErrorHandling::StdDev:
        err = Vector3dAccumulator::getStdDev();
        break;
      case ErrorHandling::Variance:
        err = Vector3dAccumulator::getVariance();
        break;
      default:
        break;
      }

      for (int i = 0; i < Vector3dAccumulator::getAverage().size(); i++) {
        if (count == 0 && errorHandling == ErrorHandling::CI95) {
          stream << "\t";
        } else {
          if (std::isinf(err[i]) || std::isnan(err[i])) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "%s",
                     errorMessage.c_str());
            painCave.isFatal = 1;
            simError();
          } else {
            if (count == 0)
              stream << "\t";
            else
              stream << "\t" << err[i];
          }
        }
      }
    }

    std::type_index getType() const override { return typeid(Vector3d); }

    std::size_t getCount() const override {
      return Vector3dAccumulator::getCount();
    }
  };

  template<>
  class AccumulatorView<PotVecAccumulator> :
      public BaseAccumulator,
      private PotVecAccumulator {
  public:
    void add(const potVec& val) override { PotVecAccumulator::add(val); }

    // Other overrides for invalid entries
    void add(RealType val) override { accumulatorFunctionCallMismatch(); }
    void add(const std::vector<RealType>&) override {
      accumulatorFunctionCallMismatch();
    }
    void add(const Vector3d&) override { accumulatorFunctionCallMismatch(); }
    void add(const Mat3x3d&) override { accumulatorFunctionCallMismatch(); }

    void writeData(std::ostream& stream, const std::string& errorMessage,
                   DataHandling dataHandling) const override {
      potVec dat;
      std::size_t count = PotVecAccumulator::getCount();

      switch (dataHandling) {
      case DataHandling::Average:
        dat = PotVecAccumulator::getAverage();
        break;
      case DataHandling::Last:
        dat = PotVecAccumulator::getLastValue();
        break;
      case DataHandling::Total:
        dat = PotVecAccumulator::getTotal();
        break;
      case DataHandling::Max:
      case DataHandling::Min:
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "%s",
                 "Max and Min functions are not defined for a "
                 "std::vector Accumulator.");
        painCave.isFatal = 1;
        simError();
        break;
      default:
        break;
      }

      for (int i = 0; i < PotVecAccumulator::getAverage().size(); i++) {
        if (std::isinf(dat[i]) || std::isnan(dat[i])) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "%s",
                   errorMessage.c_str());
          painCave.isFatal = 1;
          simError();
        } else {
          if (count == 0)
            stream << "\t";
          else
            stream << "\t" << dat[i];
        }
      }
    }

    void writeErrorBars(std::ostream& stream, const std::string& errorMessage,
                        ErrorHandling errorHandling) const override {
      potVec err;
      std::size_t count = PotVecAccumulator::getCount();

      switch (errorHandling) {
      case ErrorHandling::CI95:
        err = PotVecAccumulator::get95percentConfidenceInterval();
        break;
      case ErrorHandling::StdDev:
        err = PotVecAccumulator::getStdDev();
        break;
      case ErrorHandling::Variance:
        err = PotVecAccumulator::getVariance();
        break;
      default:
        break;
      }

      for (int i = 0; i < PotVecAccumulator::getAverage().size(); i++) {
        if (count == 0 && errorHandling == ErrorHandling::CI95) {
          stream << "\t";
        } else {
          if (std::isinf(err[i]) || std::isnan(err[i])) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "%s",
                     errorMessage.c_str());
            painCave.isFatal = 1;
            simError();
          } else {
            if (count == 0)
              stream << "\t";
            else
              stream << "\t" << err[i];
          }
        }
      }
    }

    std::type_index getType() const override { return typeid(potVec); }

    std::size_t getCount() const override {
      return PotVecAccumulator::getCount();
    }
  };

  template<>
  class AccumulatorView<Mat3x3dAccumulator> :
      public BaseAccumulator,
      private Mat3x3dAccumulator {
  public:
    void add(const Mat3x3d& val) override { Mat3x3dAccumulator::add(val); }

    // Other overrides for invalid entries
    void add(RealType val) override { accumulatorFunctionCallMismatch(); }
    void add(const std::vector<RealType>&) override {
      accumulatorFunctionCallMismatch();
    }
    void add(const Vector3d&) override { accumulatorFunctionCallMismatch(); }
    void add(const potVec&) override { accumulatorFunctionCallMismatch(); }

    void writeData(std::ostream& stream, const std::string& errorMessage,
                   DataHandling dataHandling) const override {
      Mat3x3d dat;
      std::size_t count = Mat3x3dAccumulator::getCount();

      switch (dataHandling) {
      case DataHandling::Average:
        dat = Mat3x3dAccumulator::getAverage();
        break;
      case DataHandling::Last:
        dat = Mat3x3dAccumulator::getLastValue();
        break;
      case DataHandling::Total:
        dat = Mat3x3dAccumulator::getTotal();
        break;
      case DataHandling::Max:
      case DataHandling::Min:
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "%s",
                 "Max and Min functions are not defined for a "
                 "std::vector Accumulator.");
        painCave.isFatal = 1;
        simError();
        break;
      default:
        break;
      }

      for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
          if (std::isinf(dat(i, j)) || std::isnan(dat(i, j))) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "%s",
                     errorMessage.c_str());
            painCave.isFatal = 1;
            simError();
          } else {
            if (count == 0)
              stream << "\t";
            else
              stream << "\t" << dat(i, j);
          }
        }
      }
    }

    void writeErrorBars(std::ostream& stream, const std::string& errorMessage,
                        ErrorHandling errorHandling) const override {
      Mat3x3d err;
      std::size_t count = Mat3x3dAccumulator::getCount();

      switch (errorHandling) {
      case ErrorHandling::CI95:
        err = Mat3x3dAccumulator::get95percentConfidenceInterval();
        break;
      case ErrorHandling::StdDev:
        err = Mat3x3dAccumulator::getStdDev();
        break;
      case ErrorHandling::Variance:
        err = Mat3x3dAccumulator::getVariance();
        break;
      default:
        break;
      }

      for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
          if (count == 0 && errorHandling == ErrorHandling::CI95) {
            stream << "\t";
          } else {
            if (std::isinf(err(i, j)) || std::isnan(err(i, j))) {
              snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "%s",
                       errorMessage.c_str());
              painCave.isFatal = 1;
              simError();
            } else {
              if (count == 0)
                stream << "\t";
              else
                stream << "\t" << err(i, j);
            }
          }
        }
      }
    }

    std::type_index getType() const override { return typeid(Mat3x3d); }

    std::size_t getCount() const override {
      return Mat3x3dAccumulator::getCount();
    }
  };
}  // namespace OpenMD::Utils

#endif
