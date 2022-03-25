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
 * @file GenericData.hpp
 * @brief
 * @author tlin
 * @date 09/20/2004
 * @version 1.0
 */

#ifndef UTIL_GENERICDATA_HPP
#define UTIL_GENERICDATA_HPP

#include <config.h>

#include <list>
#include <string>
#include <vector>

namespace OpenMD {

  /**
   * @ class GenericData GenericData.hpp "utils/GenericData.hpp"
   * @brief Base class for generic data which is associated with an id
   */
  class GenericData {
  public:
    GenericData() : id_("UndefinedGenericData") {}

    GenericData(const std::string& id) { setID(id); }

    /** virtual destructor */
    virtual ~GenericData() {}

    /**
     *  Returns the id of this generic data
     *
     * @return the id of this generic data
     *
     * @see #setID
     */
    const std::string getID() const { return id_; }

    /**
     *  Sets the id of this generic data
     *
     * @param id the id to be set
     *
     * @see #getID
     */
    void setID(const std::string& id) { id_ = id; }

  private:
    GenericData(const GenericData&);
    GenericData& operator=(GenericData&);
    std::string id_;
  };

  /**
   * @class SimpleTypeData
   * @brief SimpleTypeData class is a POD  repository class
   * @warning ElemDataType must be copy constructible, and copy assignable
   */
  template<typename ElemDataType>
  class SimpleTypeData : public GenericData {
  public:
    SimpleTypeData() : GenericData(), data_(ElemDataType()) {}
    SimpleTypeData(const std::string& id) :
        GenericData(id), data_(ElemDataType()) {}
    SimpleTypeData(const std::string& id, const ElemDataType& data) :
        GenericData(id), data_(data) {}
    template<typename T>
    SimpleTypeData(const SimpleTypeData<T>& s) {
      data_ = s.getData();
    }

    SimpleTypeData<ElemDataType>& operator=(
        const SimpleTypeData<ElemDataType>& s) {
      if (this == &s) return *this;

      data_ = s.getData();
      return *this;
    }

    template<typename T>
    SimpleTypeData<ElemDataType>& operator=(const SimpleTypeData<T>& s) {
      data_ = s.getData();
      return *this;
    }

    /** Returns POD data */
    const ElemDataType& getData() const { return data_; }
    ElemDataType& getData() { return data_; }
    /**
     * Sets POD data
     * @param data POD data to be set
     */
    void setData(const ElemDataType& data) { data_ = data; }

  private:
    ElemDataType data_;
  };

  /** BoolGenericData is a generic data type contains a bool variable */
  using BoolGenericData = SimpleTypeData<bool>;

  /** IntGenericData is a generic data type contains an integer variable */
  using IntGenericData = SimpleTypeData<int>;

  /** FloatGenericData is a generic data type contains a float variable */
  using FloatGenericData = SimpleTypeData<float>;

  /** DoubleGenericData is a generic data type contains a RealType variable */
  using DoubleGenericData = SimpleTypeData<RealType>;

  /**
   * @typedef StringGenericData
   * A generic data type contains a  std::string variable
   *
   * @code
   *   StringGenericData* s = new StringGenericData("MyStringGenericData");
   *   PropertyMap propMap;
   *   GenericData* gdata;
   *
   *   s->setData("OpenMD");
   *   propMap->addProperty(s);
   *
   *   gdata = propMap->getPropertyByName("MyStringGenericData");
   *   if (gdata != NULL){
   *     s = dynamic_cast<StringGenericData*>(gdata);
   *     if (s != NULL)
   *       std::cout << s->getData() << std::endl;
   *   }
   *
   * @endcode
   */
  using StringGenericData = SimpleTypeData<std::string>;

  /**
   * @class STLContainerTypeData
   * @brief STL container type generic data which is associated with an id
   *
   * \tparam ContainerType
   * \tparam ElemDataType
   */
  template<typename ElemDataType>
  class VectorTypeData : public GenericData {
  public:
    using SelfType = VectorTypeData<ElemDataType>;

    VectorTypeData(const std::string& id) : GenericData(id) {}

    VectorTypeData(const SelfType& s) : data_(s) {}

    SelfType& operator=(const SelfType& s) {
      if (this == &s) return *this;

      this->data_ = s.data_;
      return *this;
    }

  private:
    std::vector<ElemDataType> data_;
  };

  /**
   * @typedef IntVectorGenericData
   * A generic data type contains a  std::vector<int> variable.
   */
  using IntVectorGenericData = VectorTypeData<int>;

  /**
   * @typedef IntVectorGenericData
   * A generic data type contains a  std::vector<float> variable.
   */
  using FloatVectorGenericData = VectorTypeData<float>;

  /**
   * @typedef IntVectorGenericData
   * A generic data type contains a  std::vector<RealType> variable.
   */
  using DoubleVectorGenericData = VectorTypeData<RealType>;

  /**
   * @typedef StringVectorGenericData
   *  A generic data type contains a  std::vector<string> variable.
   *
   * @code
   *  StringVectorGenericData* sv = new
   * StringVectorGenericData("MyStringVector"); GenericData* gdata; PropertyMap
   * propMap; std::vector<std::string>::iterator iter;
   *
   *  sv->push_back("Hello World");
   *  sv->push_back("OpenMD");
   *
   *  propMap.addProperty(sv);
   *
   *  gdata = propMap.getPropertyByName("MyStringVector");
   *
   *  if (gdata != NULL){
   *
   *    sv = dynamic_cast<StringVectorGenericData*>(gdata);
   *
   *    if (sv != NULL){
   *      for (iter = sv->begin(); iter != sv->end(); ++ iter)
   *        std::cout << *iter << std::endl;
   *    }
   *  }
   * @endcode
   */
  using StringVectorGenericData = VectorTypeData<std::string>;
}  // namespace OpenMD

#endif  // UTIL _GENERICDATA_HPP
