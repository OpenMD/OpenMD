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

#include "applications/hydrodynamics/HydrodynamicsModel.hpp" 
#include "math/LU.hpp"
#include "math/DynamicRectMatrix.hpp"
#include "math/SquareMatrix3.hpp"
#include "utils/OOPSEConstant.hpp"
namespace oopse {
/**
 * Reference:
 * Beatriz Carrasco and Jose Gracia de la Torre, Hydrodynamic Properties of Rigid Particles:
 * Comparison of Different Modeling and Computational Procedures. 
 * Biophysical Journal, 75(6), 3044, 1999
 */

HydrodynamicsModel::HydrodynamicsModel(StuntDouble* sd, const DynamicProperty& extraParams) : sd_(sd){
    DynamicProperty::const_iterator iter;

    iter = extraParams.find("Viscosity");
    if (iter != extraParams.end()) {
        boost::any param = iter->second;
        viscosity_ = boost::any_cast<double>(param);
    }else {
        std::cout << "HydrodynamicsModel Error\n" ;
    }

    iter = extraParams.find("Temperature");
    if (iter != extraParams.end()) {
        boost::any param = iter->second;
        temperature_ = boost::any_cast<double>(param);
    }else {
        std::cout << "HydrodynamicsModel Error\n" ;
    }    
}
 
bool HydrodynamicsModel::calcHydrodyanmicsProps() {
    if (!createBeads(beads_)) {
        std::cout << "can not create beads" << std::endl;
        return false;
    }
    
    int nbeads = beads_.size();
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
                Vector3d Rij = beads_[i].pos - beads_[j].pos;
                double rij = Rij.length();
                double rij2 = rij * rij;
                double sumSigma2OverRij2 = ((beads_[i].radius*beads_[i].radius) + (beads_[i].radius*beads_[i].radius)) / rij2;                
                Mat3x3d tmpMat;
                tmpMat = outProduct(Rij, Rij) / rij2;
                double constant = 8.0 * NumericConstant::PI * viscosity_ * rij;
                Tij = ((1.0 + sumSigma2OverRij2/3.0) * I + (1.0 - sumSigma2OverRij2) * tmpMat ) / constant;
            }else {
                double constant = 1.0 / (6.0 * NumericConstant::PI * viscosity_ * beads_[i].radius);
                Tij(0, 0) = constant;
                Tij(1, 1) = constant;
                Tij(2, 2) = constant;
            }
            B.setSubMatrix(i*3, j*3, Tij);
            std::cout << Tij << std::endl;
        }
    }

    std::cout << "B=\n"
                  << B << std::endl;
    //invert B Matrix
    invertMatrix(B, C);

    std::cout << "C=\n"
                  << C << std::endl;

    //prepare U Matrix relative to arbitrary origin O(0.0, 0.0, 0.0)
    std::vector<Mat3x3d> U;
    for (int i = 0; i < nbeads; ++i) {
        Mat3x3d currU;
        currU.setupSkewMat(beads_[i].pos);
        U.push_back(currU);
    }
    
    //calculate Xi matrix at arbitrary origin O
    Mat3x3d Xitt;
    Mat3x3d Xirr;
    Mat3x3d Xitr;

    //calculate the total volume

    double volume = 0.0;
    for (std::vector<BeadParam>::iterator iter = beads_.begin(); iter != beads_.end(); ++iter) {
        volume = 4.0/3.0 * NumericConstant::PI * pow((*iter).radius,3);
    }
        
    for (std::size_t i = 0; i < nbeads; ++i) {
        for (std::size_t j = 0; j < nbeads; ++j) {
            Mat3x3d Cij;
            C.getSubMatrix(i*3, j*3, Cij);
            
            Xitt += Cij;
            Xitr += U[i] * Cij;
            Xirr += -U[i] * Cij * U[j];            
            //Xirr += -U[i] * Cij * U[j] + (0.166*6 * viscosity_ * volume) * I;            
        }
    }

    //invert Xi to get Diffusion Tensor at arbitrary origin O
    RectMatrix<double, 6, 6> Xi;    
    RectMatrix<double, 6, 6> Do;
    Xi.setSubMatrix(0, 0, Xitt);
    Xi.setSubMatrix(0, 3, Xitr.transpose());
    Xi.setSubMatrix(3, 0, Xitr);
    Xi.setSubMatrix(3, 3, Xirr);
    //invertMatrix(Xi, Do);
    double kt = OOPSEConstant::kB * temperature_ * 1.66E-2;
    //Do *= kt;    


    Mat3x3d Dott; //translational diffusion tensor at arbitrary origin O
    Mat3x3d Dorr; //rotational diffusion tensor at arbitrary origin O
    Mat3x3d Dotr; //translation-rotation couplingl diffusion tensor at arbitrary origin O

    const static Mat3x3d zeroMat(0.0);
    
    Mat3x3d XittInv(0.0);
    XittInv = Xitt.inverse();
    
    //Xirr may not be inverted,if it one of the diagonal element is zero, for example
    //( a11 a12 0)
    //( a21 a22 0)
    //( 0    0    0)
    Mat3x3d XirrInv;
    XirrInv = Xirr.inverse();

    Mat3x3d tmp;
    Mat3x3d tmpInv;
    tmp = Xitt - Xitr.transpose() * XirrInv * Xitr;
    tmpInv = tmp.inverse();

    Dott = kt * tmpInv;
    Dotr = -kt*XirrInv * Xitr * tmpInv* 1.0E8;

    tmp = Xirr - Xitr * XittInv * Xitr.transpose();    
    tmpInv = tmp.inverse();
    
    Dorr = kt * tmpInv*1.0E16;

    //Do.getSubMatrix(0, 0 , Dott);
    //Do.getSubMatrix(3, 0, Dotr);
    //Do.getSubMatrix(3, 3, Dorr);

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
    
    props_.diffCenter = rod;
    props_.transDiff = Ddtt;
    props_.transRotDiff = Ddtr;
    props_.rotDiff = Ddrr;

    return true;    
}

void HydrodynamicsModel::writeBeads(std::ostream& os) {
    std::vector<BeadParam>::iterator iter;
    os << beads_.size() << std::endl;
    os << "Generated by Hydro" << std::endl;
    for (iter = beads_.begin(); iter != beads_.end(); ++iter) {
        os << iter->atomName << "\t" << iter->pos[0] << "\t" << iter->pos[1] << "\t" << iter->pos[2] << std::endl;
    }

}

void HydrodynamicsModel::writeDiffCenterAndDiffTensor(std::ostream& os) {
    os << "//viscosity = " << viscosity_ << std::endl;
    os << "//temperature = " << temperature_<< std::endl;
    std::vector<BeadParam>::iterator iter;
    os << sd_->getType() << "\n";

    os << "//diffusion center" << std::endl;
    os << props_.diffCenter << std::endl;

    os << "//translational diffusion tensor" << std::endl;
    os << props_.transDiff << std::endl;

    os << "//translation-rotation coupling diffusion tensor" << std::endl;
    os << props_.transRotDiff << std::endl;

    os << "//rotational diffusion tensor" << std::endl;
    os << props_.rotDiff << std::endl;
    
    /*
    os << props_.diffCenter[0] << "\t" << props_.diffCenter[1] << "\t" << props_.diffCenter[2] << "\n"

    os << props_.transDiff(0, 0) << "\t" << props_.transDiff(0, 1) << "\t" << props_.transDiff(0, 2) << "\t"
        << props_.transDiff(1, 0) << "\t" << props_.transDiff(1, 1) << "\t" << props_.transDiff(1, 2) << "\t"
        << props_.transDiff(2, 0) << "\t" << props_.transDiff(2, 1) << "\t" << props_.transDiff(2, 2) << "\n";
    
    os << props_.transRotDiff(0, 0) << "\t" << props_.transRotDiff(0, 1) << "\t" << props_.transRotDiff(0, 2) << "\t"
        << props_.transRotDiff(1, 0) << "\t" << props_.transRotDiff(1, 1) << "\t" << props_.transRotDiff(1, 2) << "\t"
        << props_.transRotDiff(2, 0) << "\t" << props_.transRotDiff(2, 1) << "\t" << props_.transRotDiff(2, 2) << "\t"

    os << props_.rotDiff(0, 0) << "\t" << props_.rotDiff(0, 1) << "\t" << props_.rotDiff(0, 2) << "\t"
        << props_.rotDiff(1, 0) << "\t" << props_.rotDiff(1, 1) << "\t" << props_.rotDiff(1, 2) << "\t"
        << props_.rotDiff(2, 0) << "\t" << props_.rotDiff(2, 1) << "\t" << props_.rotDiff(2, 2) << ";"
        << std::endl;
    */
}

}
