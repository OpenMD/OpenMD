/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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

#ifndef HYDRODYNAMICS_HYDROPROP_HPP
#define HYDRODYNAMICS_HYDROPROP_HPP
#include "math/Vector3.hpp"
#include "math/SquareMatrix3.hpp"
namespace OpenMD {

  /**
   * @class HydroProp HydroProp.hpp "hydrodynamics/HydroProp.hpp"
   * Container for information about the hydrodynamic behavior of objects
   * undergoing Langevin dynamics.
   * @note the units for the center of resistance (a location) are in Angstroms
   *
   * @note the units of Xi, the resistance tensor are:
   *    Xitt (Translation-translation) : kcal fs mol^-1 Angstroms^-2 
   *    Xirt (Rotation-translation) : kcal fs mol^-1 Angstroms^-1 radians^-1 
   *    Xitr (Translation-rotation) : kcal fs mol^-1 Angstroms^-1 radians^-1 
   *    Xirr (Rotation-rotation) : kcal fs mol^-1 radians^-2 
   *
   * @note the units of D, the diffusion tensor are:
   *    Dtt (Translation-translation) : Angstroms^2 fs^-1
   *    Drt (Rotation-translation) :    Angstroms fs^-1
   *    Dtr (Translation-rotation) :    Angstroms fs^-1
   *    Drr (Rotation-rotation) :       fs^-1
   *  
   * @note after setting the value of Xi manually, the complete() function 
   * should be called to perform the Cholesky Decomposition.
   */
  class HydroProp {

  public:
    HydroProp();
    HydroProp(Vector3d cor, Mat6x6d Xi, Mat6x6d D);
    HydroProp(const std::string frictionLine);
    ~HydroProp() { } ;
    void complete();
    void setCOR(Vector3d cor) {cor_ = cor; hasCOR = true;}
    void setXi(Mat6x6d Xi) {Xi_ = Xi; hasXi = true;}
    void setD(Mat6x6d D) {D_ = D;}
    void setName(std::string name) {name_ = name;}

    Vector3d getCOR() {return cor_;}
    Mat3x3d getXitt() {return Xitt_;}
    Mat3x3d getXitr() {return Xitr_;}
    Mat3x3d getXirt() {return Xirt_;}
    Mat3x3d getXirr() {return Xirr_;}
    Mat6x6d getS() {return S_;}
    Mat6x6d getD() {return D_;}
    Mat6x6d getXi() {return Xi_;}
    std::string getName() {return name_;}    

  private:
    
    std::string name_;
    Vector3d cor_;
    Mat6x6d Xi_;
    Mat6x6d D_;
    Mat3x3d Xitt_;
    Mat3x3d Xirt_; //Xirrt == Xirtr
    Mat3x3d Xitr_;
    Mat3x3d Xirr_;
    Mat6x6d S_;
    bool hasCOR;
    bool hasXi;

  };
}

#endif
