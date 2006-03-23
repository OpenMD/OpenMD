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

#include "applications/hydrodynamics/ApproximationModel.hpp" 
#include "math/LU.hpp"
#include "math/DynamicRectMatrix.hpp"
#include "math/SquareMatrix3.hpp"
#include "utils/OOPSEConstant.hpp"
#include "applications/hydrodynamics/Spheric.hpp"
#include "applications/hydrodynamics/Ellipsoid.hpp"
#include "applications/hydrodynamics/CompositeShape.hpp"
#include "math/LU.hpp"
namespace oopse {
/**
 * Reference:
 * Beatriz Carrasco and Jose Gracia de la Torre, Hydrodynamic Properties of Rigid Particles:
 * Comparison of Different Modeling and Computational Procedures. 
 * Biophysical Journal, 75(6), 3044, 1999
 */

ApproximationModel::ApproximationModel(StuntDouble* sd, SimInfo* info): HydrodynamicsModel(sd, info){
/*
    DynamicProperty::const_iterator iter;

    iter = extraParams.find("Viscosity");
    if (iter != extraParams.end()) {
        boost::any param = iter->second;
        viscosity = boost::any_cast<double>(param);
    }else {
        std::cout << "ApproximationModel Error\n" ;
    }

    iter = extraParams.find("Temperature");
    if (iter != extraParams.end()) {
        boost::any param = iter->second;
        temperature = boost::any_cast<double>(param);
    }else {
        std::cout << "ApproximationModel Error\n" ;
    }    
*/
}

bool ApproximationModel::calcHydroProps(Spheric* spheric, double viscosity, double temperature) {
    return internalCalcHydroProps(static_cast<Shape*>(spheric), viscosity, temperature);
}

bool ApproximationModel::calcHydroProps(Ellipsoid* ellipsoid, double viscosity, double temperature) {
    return internalCalcHydroProps(static_cast<Shape*>(ellipsoid), viscosity, temperature);
}
bool ApproximationModel::calcHydroProps(CompositeShape* compositeShape, double viscosity, double temperature) {
    return internalCalcHydroProps(static_cast<Shape*>(compositeShape), viscosity, temperature);
}

 
bool ApproximationModel::internalCalcHydroProps(Shape* shape, double viscosity, double temperature) {
    if (!createBeads(beads_)) {
        std::cout << "can not create beads" << std::endl;
        return false;
    }

    bool ret = true;
    HydroProps cr;
    HydroProps cd;
    calcHydroPropsAtCR(beads_, viscosity, temperature, cr);
    calcHydroPropsAtCD(beads_, viscosity, temperature, cd);
    setCR(cr);
    setCD(cd);
    
    return true;    
}

bool ApproximationModel::calcHydroPropsAtCR(std::vector<BeadParam>& beads, double viscosity, double temperature, HydroProps& cr) {

    int nbeads = beads.size();
    DynamicRectMatrix<double> B(3*nbeads, 3*nbeads);
    DynamicRectMatrix<double> C(3*nbeads, 3*nbeads);
    Mat3x3d I;
    I(0, 0) = 1.0;
    I(1, 1) = 1.0;
    I(2, 2) = 1.0;
    
    for (std::size_t i = 0; i < nbeads; ++i) {
        for (std::size_t j = 0; j < nbeads; ++j) {
            Mat3x3d Tij;
            if (i != j ) {
                Vector3d Rij = beads[i].pos - beads[j].pos;
                double rij = Rij.length();
                double rij2 = rij * rij;
                double sumSigma2OverRij2 = ((beads[i].radius*beads[i].radius) + (beads[j].radius*beads[j].radius)) / rij2;                
                Mat3x3d tmpMat;
                tmpMat = outProduct(Rij, Rij) / rij2;
                double constant = 8.0 * NumericConstant::PI * viscosity * rij;
                Tij = ((1.0 + sumSigma2OverRij2/3.0) * I + (1.0 - sumSigma2OverRij2) * tmpMat ) / constant;
            }else {
                double constant = 1.0 / (6.0 * NumericConstant::PI * viscosity * beads[i].radius);
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

    double volume = 0.0;
    for (std::vector<BeadParam>::iterator iter = beads.begin(); iter != beads.end(); ++iter) {
        volume += 4.0/3.0 * NumericConstant::PI * pow((*iter).radius,3);
    }
        
    for (std::size_t i = 0; i < nbeads; ++i) {
        for (std::size_t j = 0; j < nbeads; ++j) {
            Mat3x3d Cij;
            C.getSubMatrix(i*3, j*3, Cij);
            
            Xiott += Cij;
            Xiotr += U[i] * Cij;
            Xiorr += -U[i] * Cij * U[j] + (6 * viscosity * volume) * I;            
        }
    }

    const double convertConstant = 6.023; //convert poise.angstrom to amu/fs
    Xiott *= convertConstant;
    Xiotr *= convertConstant;
    Xiorr *= convertConstant;
    

    
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
    

    SquareMatrix<double,6> Xir6x6;
    SquareMatrix<double,6> Dr6x6;

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
    double kt = OOPSEConstant::kB * temperature ;
    Drtt *= kt;
    Drrt *= kt;
    Drtr *= kt;
    Drrr *= kt;
    Xirtt *= OOPSEConstant::kb * temperature;
    Xirtr *= OOPSEConstant::kb * temperature;
    Xirrr *= OOPSEConstant::kb * temperature;
    

    cr.center = ror;
    cr.Xi.setSubMatrix(0, 0, Xirtt);
    cr.Xi.setSubMatrix(0, 3, Xirtr);
    cr.Xi.setSubMatrix(3, 0, Xirtr);
    cr.Xi.setSubMatrix(3, 3, Xirrr);
    cr.D.setSubMatrix(0, 0, Drtt);
    cr.D.setSubMatrix(0, 3, Drrt);
    cr.D.setSubMatrix(3, 0, Drtr);
    cr.D.setSubMatrix(3, 3, Drrr);    
    
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

bool ApproximationModel::calcHydroPropsAtCD(std::vector<BeadParam>& beads, double viscosity, double temperature, HydroProps& cr) {

    int nbeads = beads.size();
    DynamicRectMatrix<double> B(3*nbeads, 3*nbeads);
    DynamicRectMatrix<double> C(3*nbeads, 3*nbeads);
    Mat3x3d I;
    I(0, 0) = 1.0;
    I(1, 1) = 1.0;
    I(2, 2) = 1.0;
    
    for (std::size_t i = 0; i < nbeads; ++i) {
        for (std::size_t j = 0; j < nbeads; ++j) {
            Mat3x3d Tij;
            if (i != j ) {
                Vector3d Rij = beads[i].pos - beads[j].pos;
                double rij = Rij.length();
                double rij2 = rij * rij;
                double sumSigma2OverRij2 = ((beads[i].radius*beads[i].radius) + (beads[j].radius*beads[j].radius)) / rij2;                
                Mat3x3d tmpMat;
                tmpMat = outProduct(Rij, Rij) / rij2;
                double constant = 8.0 * NumericConstant::PI * viscosity * rij;
                Tij = ((1.0 + sumSigma2OverRij2/3.0) * I + (1.0 - sumSigma2OverRij2) * tmpMat ) / constant;
            }else {
                double constant = 1.0 / (6.0 * NumericConstant::PI * viscosity * beads[i].radius);
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

    double volume = 0.0;
    for (std::vector<BeadParam>::iterator iter = beads.begin(); iter != beads.end(); ++iter) {
        volume += 4.0/3.0 * NumericConstant::PI * pow((*iter).radius,3);
    }
        
    for (std::size_t i = 0; i < nbeads; ++i) {
        for (std::size_t j = 0; j < nbeads; ++j) {
            Mat3x3d Cij;
            C.getSubMatrix(i*3, j*3, Cij);
            
            Xitt += Cij;
            Xitr += U[i] * Cij;
            Xirr += -U[i] * Cij * U[j] + (6 * viscosity * volume) * I;            
        }
    }

    const double convertConstant = 6.023; //convert poise.angstrom to amu/fs
    Xitt *= convertConstant;
    Xitr *= convertConstant;
    Xirr *= convertConstant;

    double kt = OOPSEConstant::kB * temperature;

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

    SquareMatrix<double, 6> Dd;
    Dd.setSubMatrix(0, 0, Ddtt);
    Dd.setSubMatrix(0, 3, Ddtr.transpose());
    Dd.setSubMatrix(3, 0, Ddtr);
    Dd.setSubMatrix(3, 3, Ddrr);    
    SquareMatrix<double, 6> Xid;
    Ddtt *= kt;
    Ddtr *=kt;
    Ddrr *= kt;
    invertMatrix(Dd, Xid);



    //Xidtt in units of kcal*fs*mol^-1*Ang^-2
    //Xid /= OOPSEConstant::energyConvert;
    Xid *= OOPSEConstant::kb * temperature;

    cr.center = rod;
    cr.D.setSubMatrix(0, 0, Ddtt);
    cr.D.setSubMatrix(0, 3, Ddtr);
    cr.D.setSubMatrix(3, 0, Ddtr);
    cr.D.setSubMatrix(3, 3, Ddrr);
    cr.Xi = Xid;

    std::cout << "viscosity = " << viscosity << std::endl;
    std::cout << "temperature = " << temperature << std::endl;
    std::cout << "center of diffusion :" << std::endl;
    std::cout << rod << std::endl;
    std::cout << "diffusion tensor at center of diffusion " << std::endl;
    std::cout << "translation(A^2/fs) :" << std::endl;
    std::cout << Ddtt << std::endl;
    std::cout << "translation-rotation(A^3/fs):" << std::endl;
    std::cout << Ddtr << std::endl;
    std::cout << "rotation(A^4/fs):" << std::endl;
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
    std::cout << "rotation-translation (kcal*fs*mol^-1*Ang^-3):" << std::endl;
    std::cout << Xidrt << std::endl;
    std::cout << "translation-rotation(kcal*fs*mol^-1*Ang^-3):" << std::endl;
    std::cout << Xidtr << std::endl;
    std::cout << "rotation(kcal*fs*mol^-1*Ang^-4):" << std::endl;
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
