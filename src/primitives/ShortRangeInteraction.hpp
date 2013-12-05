/*
 * Copyright (c) 2013 The University of Notre Dame. All Rights Reserved.
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
 
#ifndef PRIMITIVES_SHORTRANGEINTERACTION_HPP
#define PRIMITIVES_SHORTRANGEINTERACTION_HPP

#include <vector>

#include "visitors/BaseVisitor.hpp"
#include "utils/PropertyMap.hpp"
#include "brains/Snapshot.hpp"
#include "brains/SnapshotManager.hpp"

namespace OpenMD{

  /**
   * @class ShortRangeInteraction 
   * @brief 
   *
   * A ShortRangeInteraction holds some bookeeping data for bonded
   * interactions (e.g. Bonds, Bends, Torsions, Inversions).
   *
   */
  class ShortRangeInteraction{
  public:    

    virtual ~ShortRangeInteraction();
        
    /**
     * Returns the global index of this ShortRangeInteraction.
     * @return  the global index of this ShortRangeInteraction 
     */
    int getGlobalIndex() {
      return globalIndex_;
    }

    /**
     * Sets the global index of this ShortRangeInteraction.
     * @param index new global index to be set
     */
    void setGlobalIndex(int index) {
      globalIndex_ = index;
    }
    
    /** 
     * Returns the local index of this ShortRangeInteraction 
     * @return the local index of this ShortRangeInteraction
     */
    int getLocalIndex() {
      return localIndex_;
    }

    /**
     * Sets the local index of this ShortRangeInteraction
     * @param index new index to be set
     */        
    void setLocalIndex(int index) {
      localIndex_ = index;
    }    


    /**
     * Sets the Snapshot Manager of this ShortRangeInteraction
     */
    void setSnapshotManager(SnapshotManager* sman) {
      snapshotMan_ = sman;
    }


    /** Returns the name of this ShortRangeInteraction */
    virtual std::string getName() = 0;
        
    /** Sets the name of this ShortRangeInteraction*/
    virtual void setName(const std::string& name) {}

    /**
     * <p>
     * The purpose of the Visitor Pattern is to encapsulate an
     * operation that you want to perform on the elements of a data
     * structure. In this way, you can change the operation being
     * performed on a structure without the need of changing the
     * classes of the elements that you are operating on. Using a
     * Visitor pattern allows you to decouple the classes for the data
     * structure and the algorithms used upon them
     * </p>
     * @param v visitor
     */      
    virtual void accept(BaseVisitor* v) = 0;

    /**
     * Returns the previous value of this ShortRangeInteraction
     * @return the value of this ShortRangeInteraction
     */    
    virtual RealType getPrevValue();

    /**
     * Returns the current value of this ShortRangeInteraction
     * @return the current value of this ShortRangeInteraction
     */    
    virtual RealType getValue();

    /**
     * Returns the value of this ShortRangeInteraction in specified snapshot 
     * @return the value of this ShortRangeInteraction
     * @param snapshotNo
     */    
    virtual RealType getValue(int snapshotNo) = 0;

    virtual std::vector<Atom*> getAtoms() { return atoms_; }

    /**
     * Adds property into property map
     * @param genData GenericData to be added into PropertyMap
     */
    void addProperty(GenericData* genData);

    /**
     * Removes property from PropertyMap by name
     * @param propName the name of property to be removed
     */
    void removeProperty(const std::string& propName);

    /**
     * clear all of the properties
     */
    void clearProperties();

    /**
     * Returns all names of properties
     * @return all names of properties
     */
    std::vector<std::string> getPropertyNames();

    /**
     * Returns all of the properties in PropertyMap
     * @return all of the properties in PropertyMap
     */      
    std::vector<GenericData*> getProperties();

    /**
     * Returns property 
     * @param propName name of property
     * @return a pointer point to property with propName. If no property named propName
     * exists, return NULL
     */      
    GenericData* getPropertyByName(const std::string& propName);

  protected:
        
    ShortRangeInteraction();         
    ShortRangeInteraction(const ShortRangeInteraction& sri);
    ShortRangeInteraction& operator=(const ShortRangeInteraction& sri);

    SnapshotManager* snapshotMan_;
    std::vector<Atom*> atoms_;

    int globalIndex_;
    int localIndex_;
        
  private:
    
    PropertyMap properties_;
  };

}//end namespace OpenMD
#endif //PRIMITIVES_STUNTDOUBLE_HPP
