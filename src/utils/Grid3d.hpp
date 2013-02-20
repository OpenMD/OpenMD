/*
 * Copyright (c) 2011 The University of Notre Dame. All Rights Reserved.
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
#ifndef UTILS_GRID3D_HPP
#define UTILS_GRID3D_HPP

using namespace std;
namespace OpenMD {
  /**
   * @class Grid3d
   * A generic 3d grid class
   */
  template<class Elem>
  class Grid3D {
  public:
    Grid3D(unsigned int dim1, unsigned int dim2, unsigned int dim3) : dim1_(dim1), dim2_(dim2), dim3_(dim3) {
      data_.resize(dim1_*dim2_*dim3_);
    }
    Elem& operator ()(unsigned int i, unsigned int j, unsigned int k) {
      int index = isValidGrid(i, j , k);
      assert(index != -1);
      return data_[index];
    }
    
    const Elem& operator () (unsigned int i, unsigned int j, unsigned int k) const {
      int index = isValidGrid(i, j , k);
      assert(index != -1);
      return data_[index];
    }
    
    vector<Elem> getAllNeighbors(unsigned int i, unsigned int j, unsigned int k) {
      vector<Elem> result;
      int index;
      index = isValidGrid(i-1, j, k);
      if (index != -1)
        result.push_back(data_[index]);
      
      index = isValidGrid(i+1, j, k);
      if (index != -1)
        result.push_back(data_[index]);
      
      index = isValidGrid(i, j-1, k);
      if (index != -1)
        result.push_back(data_[index]);
      
      index = isValidGrid(i, j+1, k);
      if (index != -1)
        result.push_back(data_[index]);
      
      index = isValidGrid(i, j, k-1);
      if (index != -1)
        result.push_back(data_[index]);
      
      index = isValidGrid(i, j, k+1);
      if (index != -1)
        result.push_back(data_[index]);
      
      return result;
    }

  private:
    
    int isValidGrid(unsigned int i, unsigned int j, unsigned int k) const {
      unsigned int index = i * dim2_*dim3_ + j * dim3_ + k;
      return index < data_.size() ? int(index) : -1;
    };
    
    unsigned int dim1_;
    unsigned int dim2_;
    unsigned int dim3_;
    vector<Elem> data_;
    
  };
}
#endif
