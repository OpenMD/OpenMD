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

#ifndef HYDRODYNAMICS_HYDROPROP_HPP
#define HYDRODYNAMICS_HYDROPROP_HPP
#include "math/Vector3.hpp"
#include "math/SquareMatrix3.hpp"
namespace oopse {

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
