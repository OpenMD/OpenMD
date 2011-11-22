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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include "applications/hydrodynamics/ApproximationModel.hpp" 
#include "math/LU.hpp"
#include "math/DynamicRectMatrix.hpp"
#include "math/SquareMatrix3.hpp"
#include "utils/PhysicalConstants.hpp"
#include "hydrodynamics/Sphere.hpp"
#include "hydrodynamics/Ellipsoid.hpp"
#include "applications/hydrodynamics/CompositeShape.hpp"
#include "math/LU.hpp"
#include "utils/simError.h"
namespace OpenMD {
/**
 * Reference:
 * Beatriz Carrasco and Jose Gracia de la Torre, Hydrodynamic Properties of Rigid Particles:
 * Comparison of Different Modeling and Computational Procedures. 
 * Biophysical Journal, 75(6), 3044, 1999
 */

  ApproximationModel::ApproximationModel(StuntDouble* sd, SimInfo* info): HydrodynamicsModel(sd, info){    
  }
  
  void ApproximationModel::init() {
    if (!createBeads(beads_)) {
      sprintf(painCave.errMsg, "ApproximationModel::init() : Can not create beads\n");
      painCave.isFatal = 1;
      simError();        
    }
    
  }
  
  bool ApproximationModel::calcHydroProps(Shape* shape, RealType viscosity, RealType temperature) {
    
    bool ret = true;
    HydroProp* cr = new HydroProp();
    HydroProp* cd = new HydroProp();
    calcHydroPropsAtCR(beads_, viscosity, temperature, cr);
    calcHydroPropsAtCD(beads_, viscosity, temperature, cd);
    setCR(cr);
    setCD(cd);
    return true;    
  }
  
  bool ApproximationModel::calcHydroPropsAtCR(std::vector<BeadParam>& beads, RealType viscosity, RealType temperature, HydroProp* cr) {
    
    int nbeads = beads.size();
    DynamicRectMatrix<RealType> B(3*nbeads, 3*nbeads);
    DynamicRectMatrix<RealType> C(3*nbeads, 3*nbeads);
    Mat3x3d I;
    I(0, 0) = 1.0;
    I(1, 1) = 1.0;
    I(2, 2) = 1.0;
    
    for (std::size_t i = 0; i < nbeads; ++i) {
      for (std::size_t j = 0; j < nbeads; ++j) {
        Mat3x3d Tij;
            if (i != j ) {
              Vector3d Rij = beads[i].pos - beads[j].pos;
              RealType rij = Rij.length();
              RealType rij2 = rij * rij;
              RealType sumSigma2OverRij2 = ((beads[i].radius*beads[i].radius) + (beads[j].radius*beads[j].radius)) / rij2;                
              Mat3x3d tmpMat;
              tmpMat = outProduct(Rij, Rij) / rij2;
              RealType constant = 8.0 * NumericConstant::PI * viscosity * rij;
	      RealType tmp1 = 1.0 + sumSigma2OverRij2/3.0;
	      RealType tmp2 = 1.0 - sumSigma2OverRij2;
              Tij = (tmp1 * I + tmp2 * tmpMat ) / constant;
            }else {
              RealType constant = 1.0 / (6.0 * NumericConstant::PI * viscosity * beads[i].radius);
              Tij(0, 0) = constant;
              Tij(1, 1) = constant;
              Tij(2, 2) = constant;
            }
            B.setSubMatrix(i*3, j*3, Tij);
      }
    }
    
    //invert B Matrix
    invertMatrix(B, C);
    
    //prepare U Matrix relative to arbitrary origin O(0.0, 0.0, 0.0)
    std::vector<Mat3x3d> U;
    for (int i = 0; i < nbeads; ++i) {
      Mat3x3d currU;
      currU.setupSkewMat(beads[i].pos);
      U.push_back(currU);
    }
    
    //calculate Xi matrix at arbitrary origin O
    Mat3x3d Xiott;
    Mat3x3d Xiorr;
    Mat3x3d Xiotr;
    
    //calculate the total volume
    
    RealType volume = 0.0;
    for (std::vector<BeadParam>::iterator iter = beads.begin(); iter != beads.end(); ++iter) {
      volume += 4.0/3.0 * NumericConstant::PI * pow((*iter).radius,3);
    }
    
    for (std::size_t i = 0; i < nbeads; ++i) {
      for (std::size_t j = 0; j < nbeads; ++j) {
        Mat3x3d Cij;
        C.getSubMatrix(i*3, j*3, Cij);
        
        Xiott += Cij;
        Xiotr += U[i] * Cij;
	// uncorrected here.  Volume correction is added after we assemble Xiorr
        Xiorr += -U[i] * Cij * U[j]; 
      }
    }

    // add the volume correction
    Xiorr += (6.0 * viscosity * volume) * I;    
    
    Xiott *= PhysicalConstants::viscoConvert;
    Xiotr *= PhysicalConstants::viscoConvert;
    Xiorr *= PhysicalConstants::viscoConvert;
    
    Mat3x3d tmp;
    Mat3x3d tmpInv;
    Vector3d tmpVec;
    tmp(0, 0) = Xiott(1, 1) + Xiott(2, 2);
    tmp(0, 1) = - Xiott(0, 1);
    tmp(0, 2) = -Xiott(0, 2);
    tmp(1, 0) = -Xiott(0, 1);
    tmp(1, 1) = Xiott(0, 0)  + Xiott(2, 2);
    tmp(1, 2) = -Xiott(1, 2);
    tmp(2, 0) = -Xiott(0, 2);
    tmp(2, 1) = -Xiott(1, 2);
    tmp(2, 2) = Xiott(1, 1) + Xiott(0, 0);
    tmpVec[0] = Xiotr(2, 1) - Xiotr(1, 2);
    tmpVec[1] = Xiotr(0, 2) - Xiotr(2, 0);
    tmpVec[2] = Xiotr(1, 0) - Xiotr(0, 1);
    tmpInv = tmp.inverse();    
    Vector3d ror = tmpInv * tmpVec; //center of resistance
    Mat3x3d Uor;
    Uor.setupSkewMat(ror);
    
    Mat3x3d Xirtt;
    Mat3x3d Xirrr;
    Mat3x3d Xirtr;

    Xirtt = Xiott;
    Xirtr = (Xiotr - Uor * Xiott);
    Xirrr = Xiorr - Uor * Xiott * Uor + Xiotr * Uor - Uor * Xiotr.transpose();
    

    SquareMatrix<RealType,6> Xir6x6;
    SquareMatrix<RealType,6> Dr6x6;

    Xir6x6.setSubMatrix(0, 0, Xirtt);
    Xir6x6.setSubMatrix(0, 3, Xirtr.transpose());
    Xir6x6.setSubMatrix(3, 0, Xirtr);
    Xir6x6.setSubMatrix(3, 3, Xirrr);

    invertMatrix(Xir6x6, Dr6x6);
    Mat3x3d Drtt;
    Mat3x3d Drtr;
    Mat3x3d Drrt;
    Mat3x3d Drrr;
    Dr6x6.getSubMatrix(0, 0, Drtt);
    Dr6x6.getSubMatrix(0, 3, Drrt);
    Dr6x6.getSubMatrix(3, 0, Drtr);
    Dr6x6.getSubMatrix(3, 3, Drrr);
    RealType kt = PhysicalConstants::kb * temperature ; // in kcal mol^-1
    Drtt *= kt;
    Drrt *= kt;
    Drtr *= kt;
    Drrr *= kt;
    //Xirtt *= PhysicalConstants::kb * temperature;
    //Xirtr *= PhysicalConstants::kb * temperature;
    //Xirrr *= PhysicalConstants::kb * temperature;
    
    Mat6x6d Xi, D;

    cr->setCOR(ror);

    Xi.setSubMatrix(0, 0, Xirtt);
    Xi.setSubMatrix(0, 3, Xirtr);
    Xi.setSubMatrix(3, 0, Xirtr);
    Xi.setSubMatrix(3, 3, Xirrr);

    cr->setXi(Xi);

    D.setSubMatrix(0, 0, Drtt);
    D.setSubMatrix(0, 3, Drrt);
    D.setSubMatrix(3, 0, Drtr);
    D.setSubMatrix(3, 3, Drrr);    

    cr->setD(D);
    
    std::cout << "-----------------------------------------\n";
    std::cout << "center of resistance :" << std::endl;
    std::cout << ror << std::endl;
    std::cout << "resistant tensor at center of resistance" << std::endl;
    std::cout << "translation:" << std::endl;
    std::cout << Xirtt << std::endl;
    std::cout << "translation-rotation:" << std::endl;
    std::cout << Xirtr << std::endl;
    std::cout << "rotation:" << std::endl;
    std::cout << Xirrr << std::endl; 
    std::cout << "diffusion tensor at center of resistance" << std::endl;
    std::cout << "translation:" << std::endl;
    std::cout << Drtt << std::endl;
    std::cout << "rotation-translation:" << std::endl;
    std::cout << Drrt << std::endl;
    std::cout << "translation-rotation:" << std::endl;
    std::cout << Drtr << std::endl;
    std::cout << "rotation:" << std::endl;
    std::cout << Drrr << std::endl;    
    std::cout << "-----------------------------------------\n";

    return true;
}
  
  bool ApproximationModel::calcHydroPropsAtCD(std::vector<BeadParam>& beads, RealType viscosity, RealType temperature, HydroProp* cd) {
    
    int nbeads = beads.size();
    DynamicRectMatrix<RealType> B(3*nbeads, 3*nbeads);
    DynamicRectMatrix<RealType> C(3*nbeads, 3*nbeads);
    Mat3x3d I;
    I(0, 0) = 1.0;
    I(1, 1) = 1.0;
    I(2, 2) = 1.0;
    
    for (std::size_t i = 0; i < nbeads; ++i) {
      for (std::size_t j = 0; j < nbeads; ++j) {
        Mat3x3d Tij;
        if (i != j ) {
          Vector3d Rij = beads[i].pos - beads[j].pos;
          RealType rij = Rij.length();
          RealType rij2 = rij * rij;
          RealType sumSigma2OverRij2 = ((beads[i].radius*beads[i].radius) + (beads[j].radius*beads[j].radius)) / rij2;                
          Mat3x3d tmpMat;
          tmpMat = outProduct(Rij, Rij) / rij2;
          RealType constant = 8.0 * NumericConstant::PI * viscosity * rij;
          RealType tmp1 = 1.0 + sumSigma2OverRij2/3.0;
	  RealType tmp2 = 1.0 - sumSigma2OverRij2;
	  Tij = (tmp1 * I + tmp2 * tmpMat ) / constant;
        }else {
          RealType constant = 1.0 / (6.0 * NumericConstant::PI * viscosity * beads[i].radius);
          Tij(0, 0) = constant;
          Tij(1, 1) = constant;
          Tij(2, 2) = constant;
        }
        B.setSubMatrix(i*3, j*3, Tij);
      }
    }
    
    //invert B Matrix
    invertMatrix(B, C);
    
    //prepare U Matrix relative to arbitrary origin O(0.0, 0.0, 0.0)
    std::vector<Mat3x3d> U;
    for (int i = 0; i < nbeads; ++i) {
      Mat3x3d currU;
      currU.setupSkewMat(beads[i].pos);
      U.push_back(currU);
    }
    
    //calculate Xi matrix at arbitrary origin O
    Mat3x3d Xitt;
    Mat3x3d Xirr;
    Mat3x3d Xitr;

    //calculate the total volume

    RealType volume = 0.0;
    for (std::vector<BeadParam>::iterator iter = beads.begin(); iter != beads.end(); ++iter) {
      volume += 4.0/3.0 * NumericConstant::PI * pow((*iter).radius,3);
    }
    
    for (std::size_t i = 0; i < nbeads; ++i) {
      for (std::size_t j = 0; j < nbeads; ++j) {
        Mat3x3d Cij;
        C.getSubMatrix(i*3, j*3, Cij);
            
        Xitt += Cij;
        Xitr += U[i] * Cij;
	// uncorrected here.  Volume correction is added after we assemble Xiorr
        Xirr += -U[i] * Cij * U[j];
      }
    }
    // add the volume correction here:
    Xirr += (6.0 * viscosity * volume) * I;    
    
    Xitt *= PhysicalConstants::viscoConvert;
    Xitr *= PhysicalConstants::viscoConvert;
    Xirr *= PhysicalConstants::viscoConvert;
    
    RealType kt = PhysicalConstants::kb * temperature; // in kcal mol^-1
    
    Mat3x3d Dott; //translational diffusion tensor at arbitrary origin O
    Mat3x3d Dorr; //rotational diffusion tensor at arbitrary origin O
    Mat3x3d Dotr; //translation-rotation couplingl diffusion tensor at arbitrary origin O
    
    const static Mat3x3d zeroMat(0.0);
    
    Mat3x3d XittInv(0.0);
    XittInv = Xitt.inverse();
    
    Mat3x3d XirrInv;
    XirrInv = Xirr.inverse();

    Mat3x3d tmp;
    Mat3x3d tmpInv;
    tmp = Xitt - Xitr.transpose() * XirrInv * Xitr;
    tmpInv = tmp.inverse();

    Dott = tmpInv;
    Dotr = -XirrInv * Xitr * tmpInv;
    
    tmp = Xirr - Xitr * XittInv * Xitr.transpose();    
    tmpInv = tmp.inverse();
    
    Dorr = tmpInv;

    //calculate center of diffusion
    tmp(0, 0) = Dorr(1, 1) + Dorr(2, 2);
    tmp(0, 1) = - Dorr(0, 1);
    tmp(0, 2) = -Dorr(0, 2);
    tmp(1, 0) = -Dorr(0, 1);
    tmp(1, 1) = Dorr(0, 0)  + Dorr(2, 2);
    tmp(1, 2) = -Dorr(1, 2);
    tmp(2, 0) = -Dorr(0, 2);
    tmp(2, 1) = -Dorr(1, 2);
    tmp(2, 2) = Dorr(1, 1) + Dorr(0, 0);

    Vector3d tmpVec;
    tmpVec[0] = Dotr(1, 2) - Dotr(2, 1);
    tmpVec[1] = Dotr(2, 0) - Dotr(0, 2);
    tmpVec[2] = Dotr(0, 1) - Dotr(1, 0);

    tmpInv = tmp.inverse();
    
    Vector3d rod = tmpInv * tmpVec;

    //calculate Diffusion Tensor at center of diffusion
    Mat3x3d Uod;
    Uod.setupSkewMat(rod);
    
    Mat3x3d Ddtt; //translational diffusion tensor at diffusion center
    Mat3x3d Ddtr; //rotational diffusion tensor at diffusion center
    Mat3x3d Ddrr; //translation-rotation couplingl diffusion tensor at diffusion tensor
    
    Ddtt = Dott - Uod * Dorr * Uod + Dotr.transpose() * Uod - Uod * Dotr;
    Ddrr = Dorr;
    Ddtr = Dotr + Dorr * Uod;

    SquareMatrix<RealType, 6> Dd;
    Dd.setSubMatrix(0, 0, Ddtt);
    Dd.setSubMatrix(0, 3, Ddtr.transpose());
    Dd.setSubMatrix(3, 0, Ddtr);
    Dd.setSubMatrix(3, 3, Ddrr);    
    SquareMatrix<RealType, 6> Xid;
    Ddtt *= kt;
    Ddtr *=kt;
    Ddrr *= kt;
    invertMatrix(Dd, Xid);



    //Xidtt in units of kcal*fs*mol^-1*Ang^-2
    //Xid /= PhysicalConstants::energyConvert;
    Xid *= PhysicalConstants::kb * temperature;

    Mat6x6d Xi, D;

    cd->setCOR(rod);

    cd->setXi(Xid);

    D.setSubMatrix(0, 0, Ddtt);
    D.setSubMatrix(0, 3, Ddtr);
    D.setSubMatrix(3, 0, Ddtr);
    D.setSubMatrix(3, 3, Ddrr);

    cd->setD(D);

    std::cout << "viscosity = " << viscosity << std::endl;
    std::cout << "temperature = " << temperature << std::endl;
    std::cout << "center of diffusion :" << std::endl;
    std::cout << rod << std::endl;
    std::cout << "diffusion tensor at center of diffusion " << std::endl;
    std::cout << "translation(A^2 / fs) :" << std::endl;
    std::cout << Ddtt << std::endl;
    std::cout << "translation-rotation(A / fs):" << std::endl;
    std::cout << Ddtr << std::endl;
    std::cout << "rotation(fs^-1):" << std::endl;
    std::cout << Ddrr << std::endl;

    std::cout << "resistance tensor at center of diffusion " << std::endl;
    std::cout << "translation(kcal*fs*mol^-1*Ang^-2) :" << std::endl;

    Mat3x3d Xidtt;
    Mat3x3d Xidrt;
    Mat3x3d Xidtr;
    Mat3x3d Xidrr;
    Xid.getSubMatrix(0, 0, Xidtt);
    Xid.getSubMatrix(0, 3, Xidrt);
    Xid.getSubMatrix(3, 0, Xidtr);
    Xid.getSubMatrix(3, 3, Xidrr);

    std::cout << Xidtt << std::endl;
    std::cout << "rotation-translation (kcal*fs*mol^-1*Ang^-1):" << std::endl;
    std::cout << Xidrt << std::endl;
    std::cout << "translation-rotation(kcal*fs*mol^-1*Ang^-1):" << std::endl;
    std::cout << Xidtr << std::endl;
    std::cout << "rotation(kcal*fs*mol^-1):" << std::endl;
    std::cout << Xidrr << std::endl;

    return true;
    
  }

  void ApproximationModel::writeBeads(std::ostream& os) {
    std::vector<BeadParam>::iterator iter;
    os << beads_.size() << std::endl;
    os << "Generated by Hydro" << std::endl;
    for (iter = beads_.begin(); iter != beads_.end(); ++iter) {
      os << iter->atomName << "\t" << iter->pos[0] << "\t" << iter->pos[1] << "\t" << iter->pos[2] << std::endl;
    }
    
  }    
}
