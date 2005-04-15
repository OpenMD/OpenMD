/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
 
/**
 * @file GenericData.hpp
 * @brief
 * @author tlin
 * @date 09/20/2004
 * @time 9:30am
 * @version 1.0
 */
 
#ifndef UTIL_GENERICDATA_HPP
#define UTIL_GENERICDATA_HPP

#include <list>
#include <string>
#include <vector>
namespace oopse{

  /**
   * @ class GenericData GenericData.hpp "utils/GenericData.hpp"
   * @brief Base class for generic data which is associated with an id
   */
  class GenericData{
  public:
    GenericData() :  id_("UndefinedGenericData"){}

    GenericData(const std::string& id) {  setID(id);  }

    /** virtual destructor */
    virtual ~GenericData() {}


    /**
     *  Returns the id of this generic data
     *
     * @return the id of this generic data
     *
     * @see #setID
     */
    const std::string getID() const { return id_;  }

    /**
     *  Sets the id of this generic data
     *
     * @param id the id to be set
     *
     * @see #getID
     */
    void setID(const std::string& id) { id_ = id;  }

            
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
  template<typename ElemDataType> class SimpleTypeData : public GenericData{

  public:
    SimpleTypeData() :  GenericData(), data_(ElemDataType()) {}
    SimpleTypeData(const std::string& id) : GenericData(id), data_(ElemDataType()) {}
    SimpleTypeData(const std::string&id , const ElemDataType& data) : GenericData(id), data_(data) {}
    template<typename T>
    SimpleTypeData(const SimpleTypeData<T>& s) {
      data_ = s.getData();
    }

    SimpleTypeData<ElemDataType>& operator =(const SimpleTypeData<ElemDataType>& s) {
      if (this == &s)
	return *this;
                
      data_ = s.getData();
      return *this;                
    }
            
    template<typename T>
    SimpleTypeData<ElemDataType>& operator =(const SimpleTypeData<T>& s) {                
      data_ = s.getData();
      return *this;
    }
            
    /** Returns POD data */     
    const ElemDataType& getData() const {return data_;}
    ElemDataType& getData()  {return data_;}
    /**
     * Sets POD data
     * @data POD data to be set
     */
    void setData(const ElemDataType& data) { data_ = data;  }

  private:
    ElemDataType data_;
  };

  /** BoolGenericData is a generic data type contains a bool variable */
  typedef SimpleTypeData<bool> BoolGenericData;

  /** IntGenericData is a generic data type contains an integer variable */
  typedef SimpleTypeData<int> IntGenericData;

  /** FloatGenericData is a generic data type contains a float variable */
  typedef SimpleTypeData<float> FloatGenericData;

  /** DoubleGenericData is a generic data type contains a double variable */
  typedef SimpleTypeData<double> DoubleGenericData;
  
  /**
   * @typedef StringGenericData
   * A generic data type contains a  std::string variable
   *
   * @code
   *   StringGenericData* s = new StringGenericData("MyStringGenericData");
   *   PropertyMap propMap;
   *   GenericData* gdata;
   *
   *   s->setData("OOPSE");
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
  typedef SimpleTypeData<std::string> StringGenericData;

  /**
   * @class STLContainerTypeData 
   * @brief STL container type generic data which is associated with an id
   *
   * @template ContainerType
   * @template ElemDataType
   */
  template <typename ElemDataType > 
  class VectorTypeData : public GenericData {
  public:
    typedef VectorTypeData<ElemDataType> SelfType;

    VectorTypeData(const std::string& id) 
      : GenericData(id){}
            
    VectorTypeData(const SelfType& s) : SelfType(s){}

    SelfType& operator =(const SelfType& s){
      if (this == &s)
	return *this;

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
  typedef VectorTypeData<int> IntVectorGenericData;

  /**
   * @typedef IntVectorGenericData
   * A generic data type contains a  std::vector<float> variable.
   */  
  typedef VectorTypeData<float> FloatVectorGenericData;

  /**
   * @typedef IntVectorGenericData
   * A generic data type contains a  std::vector<double> variable.
   */  
  typedef VectorTypeData<double> DoubleVectorGenericData;

  /** 
   * @typedef StringVectorGenericData
   *  A generic data type contains a  std::vector<string> variable.
   *
   * @code
   *  StringVectorGenericData* sv = new StringVectorGenericData("MyStringVector");
   *  GenericData* gdata;
   *  PropertyMap propMap;
   *  std::vector<std::string>::iterator iter;
   *  
   *  sv->push_back("Hello World");
   *  sv->push_back("OOPSE");
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
  typedef VectorTypeData<std::string> StringVectorGenericData;
  

} // namespace oopse
#endif //UTIL _GENERICDATA_HPP
