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

#ifndef UTILS_GRID3D_HPP
#define UTILS_GRID3D_HPP

#include <cassert>
#include <vector>

namespace OpenMD {

  /**
   * @class Grid3d
   * A generic 3d grid class
   */
  template<class Elem>
  class Grid3D {
  public:
    Grid3D(unsigned int dim1, unsigned int dim2, unsigned int dim3) :
        dim1_(dim1), dim2_(dim2), dim3_(dim3) {
      data_.resize(dim1_ * dim2_ * dim3_);
    }

    Elem& operator()(unsigned int i, unsigned int j, unsigned int k) {
      int index = isValidGrid(i, j, k);
      assert(index != -1);
      return data_[index];
    }

    const Elem& operator()(unsigned int i, unsigned int j,
                           unsigned int k) const {
      int index = isValidGrid(i, j, k);
      assert(index != -1);
      return data_[index];
    }

    std::vector<Elem> getAllNeighbors(unsigned int i, unsigned int j,
                                      unsigned int k) {
      std::vector<Elem> result;
      int index;
      index = isValidGrid(i - 1, j, k);
      if (index != -1) result.push_back(data_[index]);

      index = isValidGrid(i + 1, j, k);
      if (index != -1) result.push_back(data_[index]);

      index = isValidGrid(i, j - 1, k);
      if (index != -1) result.push_back(data_[index]);

      index = isValidGrid(i, j + 1, k);
      if (index != -1) result.push_back(data_[index]);

      index = isValidGrid(i, j, k - 1);
      if (index != -1) result.push_back(data_[index]);

      index = isValidGrid(i, j, k + 1);
      if (index != -1) result.push_back(data_[index]);

      return result;
    }

  private:
    int isValidGrid(unsigned int i, unsigned int j, unsigned int k) const {
      unsigned int index = i * dim2_ * dim3_ + j * dim3_ + k;
      return index < data_.size() ? int(index) : -1;
    };

    unsigned int dim1_;
    unsigned int dim2_;
    unsigned int dim3_;
    std::vector<Elem> data_;
  };
}  // namespace OpenMD

#endif
