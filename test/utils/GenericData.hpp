/*
 * Copyright (C) 2000-2004  Object Oriented Parallel Simulation Engine (OOPSE) project
 * 
 * Contact: oopse@oopse.org
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
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
     * A generic data type contains a string variable
     *
     * @code
     *   StringGenericData* s = new StringGenericData("MyStringGenericData");
     *   PropertyMap propMap;
     *   GenericData* gdata;
     *
     *   s->setData("OOPSE");
     *   propMap->addProperty(s);
     *   
     *   gdata = propMap->getProperty("MyStringGenericData");
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
    template <template<typename ELEM, typename = std::allocator<ELEM> > class ContainerType,
                     typename ElemDataType > 
    class STLContainerTypeData : public GenericData, public ContainerType<ElemDataType>{
        public:
            typedef STLContainerTypeData<ContainerType, ElemDataType> SelfType;
            typedef ContainerType<ElemDataType> STLContainerType;

            STLContainerTypeData(const std::string& id) 
                : GenericData(id),  ContainerType<ElemDataType> () {}
            
            STLContainerTypeData(const SelfType& s) : SelfType(s){}

            SelfType& operator =(const SelfType& s){
                if (this == &s)
                    return *this;

                STLContainerType::operator=(s);
                return *this;
            }
    };

    /**
     * @typedef IntVectorGenericData
     * A generic data type contains a vector<int> variable.
     */  
    typedef STLContainerTypeData<std::vector, int> IntVectorGenericData;

    /**
     * @typedef IntVectorGenericData
     * A generic data type contains a vector<float> variable.
     */  
    typedef STLContainerTypeData<std::vector, float> FloatVectorGenericData;

    /**
     * @typedef IntVectorGenericData
     * A generic data type contains a vector<double> variable.
     */  
    typedef STLContainerTypeData<std::vector, double> DoubleVectorGenericData;

    /** 
     * @typedef StringVectorGenericData
     *  A generic data type contains a vector<string> variable.
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
     *  gdata = propMap.getProperty("MyStringVector");
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
    typedef STLContainerTypeData<std::vector, std::string> StringVectorGenericData;

    /**
     * @typedef IntVectorGenericData
     * A generic data type contains a  list<vector<string> >  variable.
     */  
    typedef STLContainerTypeData<std::list, std::vector<int> > IntVectorListGenericData;
  
} // namespace oopse
#endif //UTIL _GENERICDATA_HPP
