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

#ifndef PRIMITIVES_SHORTRANGEINTERACTION_HPP
#define PRIMITIVES_SHORTRANGEINTERACTION_HPP

#include <memory>
#include <vector>

#include "brains/Snapshot.hpp"
#include "brains/SnapshotManager.hpp"
#include "utils/PropertyMap.hpp"
#include "visitors/BaseVisitor.hpp"

namespace OpenMD {

  /**
   * @class ShortRangeInteraction
   * @brief
   *
   * A ShortRangeInteraction holds some bookeeping data for bonded
   * interactions (e.g. Bonds, Bends, Torsions, Inversions).
   *
   */
  class ShortRangeInteraction {
  public:
    virtual ~ShortRangeInteraction();

    /**
     * Returns the global index of this ShortRangeInteraction.
     * @return  the global index of this ShortRangeInteraction
     */
    int getGlobalIndex() { return globalIndex_; }

    /**
     * Sets the global index of this ShortRangeInteraction.
     * @param index new global index to be set
     */
    void setGlobalIndex(int index) { globalIndex_ = index; }

    /**
     * Returns the local index of this ShortRangeInteraction
     * @return the local index of this ShortRangeInteraction
     */
    int getLocalIndex() { return localIndex_; }

    /**
     * Sets the local index of this ShortRangeInteraction
     * @param index new index to be set
     */
    void setLocalIndex(int index) { localIndex_ = index; }

    /**
     * Sets the Snapshot Manager of this ShortRangeInteraction
     */
    void setSnapshotManager(SnapshotManager* sman) { snapshotMan_ = sman; }

    /** Returns the name of this ShortRangeInteraction */
    virtual std::string getName() = 0;

    /** Sets the name of this ShortRangeInteraction*/
    virtual void setName(const std::string&) {}

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
    void addProperty(std::shared_ptr<GenericData> genData);

    /**
     * Removes property from PropertyMap by name
     * @param propName the name of property to be removed
     */
    void removeProperty(const std::string& propName);

    /**
     * Returns all names of properties
     * @return all names of properties
     */
    std::vector<std::string> getPropertyNames();

    /**
     * Returns all of the properties in PropertyMap
     * @return all of the properties in PropertyMap
     */
    std::vector<std::shared_ptr<GenericData>> getProperties();

    /**
     * Returns property
     * @param propName name of property
     * @return a pointer point to property with propName. If no property named
     * propName exists, return NULL
     */
    std::shared_ptr<GenericData> getPropertyByName(const std::string& propName);

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
}  // namespace OpenMD

#endif  // PRIMITIVES_STUNTDOUBLE_HPP
