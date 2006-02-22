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
namespace oopse {
/**
 * Reference:
 * Beatriz Carrasco and Jose Gracia de la Torre, Hydrodynamic Properties of Rigid Particles:
 * Comparison of Different Modeling and Computational Procedures. 
 * Biophysical Journal, 75(6), 3044, 1999
 */
bool HydrodynamicsModel::calcHydrodyanmicsProps(double eta) {
    if (!createBeads(beads_)) {
        std::cout << "can not create beads" << std::endl;
        return false;
    }
    
    int nbeads = beads_.size();
    DynamicRectMatrix<double> B(3*nbeads, 3*nbeads);
    DynamicRectMatrix<double> C(3*nbeads, 3*nbeads);
    Mat3x3d I;
    for (std::size_t i = 0; i < nbeads; ++i) {
        for (std::size_t j = 0; j < nbeads; ++j) {
            Mat3x3d Tij;
            if (i != j ) {
                Vector3d Rij = beads_[i].pos - beads_[j].pos;
                double rij = Rij.length();
                double rij2 = rij * rij;
                double sumSigma2OverRij2 = ((beads_[i].radius*beads_[i].radius) + (beads_[i].radius*beads_[i].radius)) / rij2;                
                Mat3x3d tmpMat;
                tmpMat = outProduct(beads_[i].pos, beads_[j].pos) / rij2;
                double constant = 8.0 * NumericConstant::PI * eta * rij;
                Tij = ((1.0 + sumSigma2OverRij2/3.0) * I + (1.0 - sumSigma2OverRij2) * tmpMat ) / constant;
            }else {
                double constant = 1.0 / (6.0 * NumericConstant::PI * eta * beads_[i].radius);
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
        currU.setupSkewMat(beads_[i].pos);
        U.push_back(currU);
    }
    
    //calculate Xi matrix at arbitrary origin O
    Mat3x3d Xitt;
    Mat3x3d Xirr;
    Mat3x3d Xitr;
        
    for (std::size_t i = 0; i < nbeads; ++i) {
        for (std::size_t j = 0; j < nbeads; ++j) {
            Mat3x3d Cij;
            C.getSubMatrix(i*3, j*3, Cij);
            
            Xitt += Cij;
            Xirr += U[i] * Cij;
            Xitr += U[i] * Cij * U[j];            
        }
    }

    //invert Xi to get Diffusion Tensor at arbitrary origin O
    RectMatrix<double, 6, 6> Xi;    
    RectMatrix<double, 6, 6> Do;
    Xi.setSubMatrix(0, 0, Xitt);
    Xi.setSubMatrix(0, 3, Xitr.transpose());
    Xi.setSubMatrix(3, 0, Xitr);
    Xi.setSubMatrix(3, 3, Xitt);
    invertMatrix(Xi, Do);

    Mat3x3d Dott; //translational diffusion tensor at arbitrary origin O
    Mat3x3d Dorr; //rotational diffusion tensor at arbitrary origin O
    Mat3x3d Dotr; //translation-rotation couplingl diffusion tensor at arbitrary origin O
    Do.getSubMatrix(0, 0 , Dott);
    Do.getSubMatrix(3, 0, Dotr);
    Do.getSubMatrix(3, 3, Dorr);

    //calculate center of diffusion
    Mat3x3d tmpMat;
    tmpMat(0, 0) = Dorr(1, 1) + Dorr(2, 2);
    tmpMat(0, 1) = - Dorr(0, 1);
    tmpMat(0, 2) = -Dorr(0, 2);
    tmpMat(1, 0) = -Dorr(0, 1);
    tmpMat(1, 1) = Dorr(0, 0)  + Dorr(2, 2);
    tmpMat(1, 2) = -Dorr(1, 2);
    tmpMat(2, 0) = -Dorr(0, 2);
    tmpMat(2, 1) = -Dorr(1, 2);
    tmpMat(2, 2) = Dorr(1, 1) + Dorr(0, 0);

    Vector3d tmpVec;
    tmpVec[0] = Dotr(1, 2) - Dotr(2, 1);
    tmpVec[1] = Dotr(2, 0) - Dotr(0, 2);
    tmpVec[2] = Dotr(0, 1) - Dotr(1, 0);
        
    Vector3d rod = tmpMat.inverse() * tmpVec;

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

}

void HydrodynamicsModel::writeDiffCenterAndDiffTensor(std::ostream& os) {

}

}
