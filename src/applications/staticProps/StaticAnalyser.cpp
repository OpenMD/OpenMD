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

#include "applications/staticProps/StaticAnalyser.hpp"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "brains/SimInfo.hpp"
#include "math/Vector3.hpp"
#include "utils/Accumulator.hpp"
#include "utils/Revision.hpp"
#include "utils/simError.h"

namespace OpenMD {

  StaticAnalyser::StaticAnalyser(SimInfo* info, const std::string& filename,
                                 unsigned int nbins) :
      info_(info),
      dumpFilename_(filename), step_(1), nBins_(nbins) {}

  void StaticAnalyser::writeOutput() {
    std::vector<OutputData*>::iterator i;
    OutputData* outputData;

    std::ofstream ofs(outputFilename_.c_str());

    if (ofs.is_open()) {
      Revision r;
      ofs << "# " << getAnalysisType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      // ofs << "# selection script1: \"" << selectionScript1_ ;
      // ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";
      ofs << "#";
      for (outputData = beginOutputData(i); outputData;
           outputData = nextOutputData(i)) {
        ofs << "\t" << outputData->title << "(" << outputData->units << ")";
        // add some extra tabs for column alignment
        if (outputData->dataType == odtVector3) ofs << "\t\t";
        if (outputData->dataType == odtArray2d) {
          ofs << "(";
          for (unsigned int j = 0; j < outputData->accumulatorArray2d[0].size();
               j++) {
            ofs << outputData->columnNames[j] << "\t";
          }
          ofs << ")\t";
        }
      }

      ofs << std::endl;

      ofs.precision(8);

      for (unsigned int j = 0; j < nBins_; j++) {
        for (outputData = beginOutputData(i); outputData;
             outputData = nextOutputData(i)) {
          writeData(ofs, outputData, j);
        }
        ofs << std::endl;
      }

      ofs << "#######################################################\n";
      ofs << "# 95% confidence intervals in those quantities follow:\n";
      ofs << "#######################################################\n";

      for (unsigned int j = 0; j < nBins_; j++) {
        ofs << "#";
        for (outputData = beginOutputData(i); outputData;
             outputData = nextOutputData(i)) {
          writeErrorBars(ofs, outputData, j);
        }
        ofs << std::endl;
      }

      ofs.flush();
      ofs.close();

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "StaticAnalyser: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }
  }

  void StaticAnalyser::writeReal(std::ostream& os, OutputData* dat,
                                 unsigned int bin) {
    assert(bin < nBins_);

    RealType r;
    std::size_t count = dat->accumulator[bin]->count();

    if (count == 0) {
      os << "\t";
    } else {
      if (dat->dataHandling == odhMax) {
        dynamic_cast<Accumulator*>(dat->accumulator[bin])->getMax(r);
      } else if (dat->dataHandling == odhTotal) {
        dynamic_cast<Accumulator*>(dat->accumulator[bin])->getTotal(r);
      } else {
        dynamic_cast<Accumulator*>(dat->accumulator[bin])->getAverage(r);
      }

      if (std::isinf(r) || std::isnan(r)) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "StaticAnalyser detected a numerical error writing:\n"
                 "\t%s for bin %u",
                 dat->title.c_str(), bin);
        painCave.isFatal = 1;
        simError();
      } else {
        os << "\t" << r;
      }
    }
  }

  void StaticAnalyser::writeVector(std::ostream& os, OutputData* dat,
                                   unsigned int bin) {
    assert(bin < nBins_);

    Vector3d v;
    std::size_t count = dat->accumulator[bin]->count();

    if (count == 0) {
      os << "\t\t\t";
    } else {
      if (dat->dataHandling == odhTotal) {
        dynamic_cast<VectorAccumulator*>(dat->accumulator[bin])->getTotal(v);
      } else {
        dynamic_cast<VectorAccumulator*>(dat->accumulator[bin])->getAverage(v);
      }

      if (std::isinf(v[0]) || std::isnan(v[0]) || std::isinf(v[1]) ||
          std::isnan(v[1]) || std::isinf(v[2]) || std::isnan(v[2])) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "StaticAnalyser detected a numerical error writing:\n"
                 "\t%s for bin %u",
                 dat->title.c_str(), bin);
        painCave.isFatal = 1;
        simError();
      } else {
        os << "\t" << v[0] << "\t" << v[1] << "\t" << v[2];
      }
    }
  }

  void StaticAnalyser::writeArray(std::ostream& os, OutputData* dat,
                                  unsigned int bin) {
    assert(bin < nBins_);

    RealType s;
    std::size_t columns = dat->accumulatorArray2d[0].size();

    for (std::size_t i = 0; i < columns; i++) {
      std::size_t count = dat->accumulatorArray2d[bin][i]->count();

      if (count == 0) {
        os << "\t";
      } else {
        if (dat->dataHandling == odhTotal) {
          dynamic_cast<Accumulator*>(dat->accumulatorArray2d[bin][i])
              ->getTotal(s);
        } else {
          dynamic_cast<Accumulator*>(dat->accumulatorArray2d[bin][i])
              ->getAverage(s);
        }

        if (std::isinf(s) || std::isnan(s)) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "StaticAnalyser detected a numerical error writing:\n"
                   "\t%s for bin %u, column %u",
                   dat->title.c_str(), bin, static_cast<unsigned int>(i));
          painCave.isFatal = 1;
          simError();
        } else {
          os << "\t" << s;
        }
      }
    }
  }

  void StaticAnalyser::writeData(std::ostream& os, OutputData* dat,
                                 unsigned int bin) {
    assert(bin < nBins_);

    if (dat->dataType == odtReal) {
      writeReal(os, dat, bin);
    } else if (dat->dataType == odtVector3) {
      writeVector(os, dat, bin);
    } else if (dat->dataType == odtArray2d) {
      writeArray(os, dat, bin);
    }
  }

  void StaticAnalyser::writeRealErrorBars(std::ostream& os, OutputData* dat,
                                          unsigned int bin) {
    assert(bin < nBins_);

    RealType r;
    std::size_t count = dat->accumulator[bin]->count();

    if (count == 0) {
      os << "\t";
    } else {
      dynamic_cast<Accumulator*>(dat->accumulator[bin])
          ->get95percentConfidenceInterval(r);

      if (std::isinf(r) || std::isnan(r)) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "StaticAnalyser detected a numerical error writing:\n"
                 "\tstandard deviation of %s for bin %u",
                 dat->title.c_str(), bin);
        painCave.isFatal = 1;
        simError();
      } else {
        os << "\t" << r;
      }
    }
  }

  void StaticAnalyser::writeVectorErrorBars(std::ostream& os, OutputData* dat,
                                            unsigned int bin) {
    assert(bin < nBins_);

    Vector3d v;
    std::size_t count = dat->accumulator[bin]->count();

    if (count == 0) {
      os << "\t\t\t";
    } else {
      dynamic_cast<VectorAccumulator*>(dat->accumulator[bin])
          ->get95percentConfidenceInterval(v);

      if (std::isinf(v[0]) || std::isnan(v[0]) || std::isinf(v[1]) ||
          std::isnan(v[1]) || std::isinf(v[2]) || std::isnan(v[2])) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "StaticAnalyser detected a numerical error writing:\n"
                 "\tstandard deviation of %s for bin %u",
                 dat->title.c_str(), bin);
        painCave.isFatal = 1;
        simError();
      } else {
        os << "\t" << v[0] << "\t" << v[1] << "\t" << v[2];
      }
    }
  }

  void StaticAnalyser::writeArrayErrorBars(std::ostream& os, OutputData* dat,
                                           unsigned int bin) {
    assert(bin < nBins_);

    RealType s;
    std::size_t columns = dat->accumulatorArray2d[0].size();

    for (std::size_t i = 0; i < columns; i++) {
      std::size_t count = dat->accumulatorArray2d[bin][i]->count();

      if (count == 0) {
        os << "\t";
      } else {
        dynamic_cast<Accumulator*>(dat->accumulatorArray2d[bin][i])
            ->get95percentConfidenceInterval(s);

        if (std::isinf(s) || std::isnan(s)) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "StaticAnalyser detected a numerical error writing:\n"
                   "\t%s std. dev. for bin %u, column %u",
                   dat->title.c_str(), bin, static_cast<unsigned int>(i));
          painCave.isFatal = 1;
          simError();
        } else {
          os << "\t" << s;
        }
      }
    }
  }

  void StaticAnalyser::writeErrorBars(std::ostream& os, OutputData* dat,
                                      unsigned int bin) {
    assert(bin < nBins_);

    if (dat->dataType == odtReal) {
      writeRealErrorBars(os, dat, bin);
    } else if (dat->dataType == odtVector3) {
      writeVectorErrorBars(os, dat, bin);
    } else if (dat->dataType == odtArray2d) {
      writeArrayErrorBars(os, dat, bin);
    }
  }

  OutputData* StaticAnalyser::beginOutputData(
      std::vector<OutputData*>::iterator& i) {
    i = data_.begin();
    return i != data_.end() ? *i : NULL;
  }

  OutputData* StaticAnalyser::nextOutputData(
      std::vector<OutputData*>::iterator& i) {
    ++i;
    return i != data_.end() ? *i : NULL;
  }
}  // namespace OpenMD
