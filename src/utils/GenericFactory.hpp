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
 * @file Vector3.hpp
 * @author Teng Lin
 * @date 09/14/2004
 * @version 1.0
 */
#ifndef UTIL_GENERICFACTORY_HPP
#define UTIL_GENERICFACTORY_HPP
#include <map>
#include <string>
namespace oopse {
template<class Product, typename ProductIdentType = std::string, typename Creator = Product* (*)() >
class GenericFactory{
    public:
        
        typedef std::map<ProductIdentType,  Creator*> CreatorMapType;
        typedef GenericFactory<Product, ProductIdentType, Creator> SelfType;
        static SelfType* getInstance() {
            if (instance_ == NULL) {
                instance_ = new GenericFactory<Product, ProductIdentType, Creator>();
            }

            return instance_;
        }
        
        //bool register( const ProductIdentType& id, Creator creator) {

            //insert method in std::map will return a pair<iterator, bool>. the second
            //element of this pair indicates whether the result of the operation
            //return creatorMap_.insert(CreatorMapType::value_type(id, creator)).second;
        //}

        bool unregister(const ProductIdentType& id) {
            
            return creatorMap_->erase(id) == 1; 

        }
        
        bool hasCreator( const ProductIdentType& id ) {
            CreatorMapType::iterator i;

            i = creatorMap_.find(id);

            if (i != creatorMap_.end()) {
                return true;
            } else {
                return false;
            }
        }

        //const std::string toString() {

        
        //}

        Product* createProduct( const ProductIdentType& id ) {
            CreatorMapType::iterator i;

            i = creatorMap_.find(id);

            if (i != creatorMap_.end()) {
                //call the function to create the product
                return (i->second)();
            } else {
                return NULL;
            }
        }
        
    private:
        GenericFactory(){}
        static SelfType* instance_;
        CreatorMapType creatorMap_;
};


#define REGISTER_CREATOR(factory, product, id) \
    product * create##product() {\
        return new product(); \
    }

}//namespace oopse
#endif //UTIL_GENERICFACTORY_HPP

