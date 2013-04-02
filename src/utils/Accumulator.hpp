/*
 * Copyright (c) 2012 The University of Notre Dame. All Rights Reserved.
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

#ifndef UTILS_ACCUMULATOR_HPP
#define UTILS_ACCUMULATOR_HPP

#include <cmath>
#include <cassert>
#include "math/Vector3.hpp"

namespace OpenMD {


  class BaseAccumulator {
  public:
    virtual void clear() = 0;
    /**
     * get the number of accumulated values
     */
    virtual size_t count()  {
      return Count_;
    }
  protected:
    size_t Count_;

  };   



  /** 
   * Basic Accumulator class for numbers. 
   */  
  class Accumulator : public BaseAccumulator {    

    typedef RealType ElementType;
    typedef RealType ResultType;

  public:
    
    Accumulator() : BaseAccumulator() {
      this->clear();
    }

    /**
     * Accumulate another value
     */
    virtual void add(ElementType const& val) {
      Count_++;
      Avg_  += (val       - Avg_ ) / Count_;
      Avg2_ += (val * val - Avg2_) / Count_;
      Val_   = val;
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
      Val_   = 0;
    }
    

    /**
     * return the most recently added value
     */
    void getLastValue(ElementType &ret)  {
      ret = Val_;
      return;
    }    

    /**
     * compute the Mean
     */
    void getAverage(ResultType &ret)  {
      assert(Count_ != 0);
      ret = Avg_;
      return;
    }

    /**
     * compute the Variance
     */
    void getVariance(ResultType &ret)  {
      assert(Count_ != 0);
      ret = (Avg2_ - Avg_  * Avg_);
      return;
    }
    
    /**
     * compute error of average value
     */
    void getStdDev(ResultType &ret)  {
      assert(Count_ != 0);
      RealType var;
      this->getVariance(var);
      ret = sqrt(var);
      return;
    }

    /**
     * return the largest value
     */
    void getMax(ElementType &ret)  {
      assert(Count_ != 0);
      ret = Max_;
      return;
    }

    /**
     * return the smallest value
     */
    void getMin(ElementType &ret)  {
      assert(Count_ != 0);
      ret = Max_;
      return;
    }

  private:
    ElementType Val_;
    ResultType Avg_;
    ResultType Avg2_;
    ElementType Min_;
    ElementType Max_;
  };

  class VectorAccumulator : public BaseAccumulator {
    
    typedef Vector3d ElementType;
    typedef Vector3d ResultType;
    
  public:
    VectorAccumulator() : BaseAccumulator() {
      this->clear();
    }

    /**
     * Accumulate another value
     */
    void add(ElementType const& val) {
      Count_++;
      RealType len(0.0);
      for (unsigned int i =0; i < 3; i++) {
        Avg_[i]  += (val[i]       - Avg_[i] ) / Count_;
        Avg2_[i] += (val[i] * val[i] - Avg2_[i]) / Count_;
        Val_[i]   = val[i];
        len += val[i]*val[i];
      }
      len = sqrt(len);
      AvgLen_  += (len       - AvgLen_ ) / Count_;
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
      Avg_ = V3Zero;
      Avg2_ = V3Zero;
      Val_ = V3Zero;
      AvgLen_   = 0;
      AvgLen2_  = 0;
    }
    
    /**
     * return the most recently added value
     */
    void getLastValue(ElementType &ret) {
      ret = Val_;
      return;
    }
    
    /**
     * compute the Mean
     */
    void getAverage(ResultType &ret) {
      assert(Count_ != 0);
      ret = Avg_;
      return;
    }
    
    /**
     * compute the Variance
     */
    void getVariance(ResultType &ret) {
      assert(Count_ != 0);
      for (unsigned int i =0; i < 3; i++) {
        ret[i] = (Avg2_[i] - Avg_[i]  * Avg_[i]);
      }
      return;
    }
    
    /**
     * compute error of average value
     */
    void getStdDev(ResultType &ret) {
      assert(Count_ != 0);
      ResultType var;
      this->getVariance(var);
      ret[0] = sqrt(var[0]);
      ret[1] = sqrt(var[1]);
      ret[2] = sqrt(var[2]);
      return;
    }

    /**
     * return the largest length
     */
    void getMaxLength(RealType &ret) {
      assert(Count_ != 0);
      ret = Max_;
      return;
    }

    /**
     * return the smallest length
     */
    void getMinLength(RealType &ret) {
      assert(Count_ != 0);
      ret = Min_;
      return;
    }

    /**
     * return the largest length
     */
    void getAverageLength(RealType &ret) {
      assert(Count_ != 0);
      ret = AvgLen_;
      return;
    }

    /**
     * compute the Variance of the length
     */
    void getLengthVariance(RealType &ret) {
      assert(Count_ != 0);      
      ret= (AvgLen2_ - AvgLen_ * AvgLen_);
      return;
    }
    
    /**
     * compute error of average value
     */
    void getLengthStdDev(RealType &ret) {
      assert(Count_ != 0);
      RealType var;
      this->getLengthVariance(var);
      ret = sqrt(var);
      return;
    }

  private:
    ResultType Val_;
    ResultType Avg_;
    ResultType Avg2_;
    RealType AvgLen_;
    RealType AvgLen2_;
    RealType Min_;
    RealType Max_;

  };

  class MatrixAccumulator : public BaseAccumulator {
    
    typedef Mat3x3d ElementType;
    typedef Mat3x3d ResultType;
    
  public:
    MatrixAccumulator() : BaseAccumulator() {
      this->clear();
    }

    /**
     * Accumulate another value
     */
    void add(ElementType const& val) {
      Count_++;
      for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {          
          Avg_(i,j)  += (val(i,j)       - Avg_(i,j) ) / Count_;
          Avg2_(i,j) += (val(i,j) * val(i,j) - Avg2_(i,j)) / Count_;
          Val_(i,j)   = val(i,j);
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
    }
    
    /**
     * return the most recently added value
     */
    void getLastValue(ElementType &ret) {
      ret = Val_;
      return;
    }
    
    /**
     * compute the Mean
     */
    void getAverage(ResultType &ret) {
      assert(Count_ != 0);
      ret = Avg_;
      return;
    }

    /**
     * compute the Variance
     */
    void getVariance(ResultType &ret) {
      assert(Count_ != 0);
      for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {          
          ret(i,j) = (Avg2_(i,j) - Avg_(i,j)  * Avg_(i,j));
        }
      }
      return;
    }
    
    /**
     * compute error of average value
     */
    void getStdDev(ResultType &ret) {
      assert(Count_ != 0);
      Mat3x3d var;
      this->getVariance(var);
      for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
          ret(i,j) = sqrt(var(i,j));  
        }
      }
      return;
    }
        
  private:
    ElementType Val_;
    ResultType Avg_;
    ResultType Avg2_;
  };


} 

#endif
