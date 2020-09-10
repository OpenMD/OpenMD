/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include "applications/hydrodynamics/ApproximationModel.hpp"
#include "math/LU.hpp"
#include "math/DynamicRectMatrix.hpp"
#include "math/SquareMatrix3.hpp"
#include "utils/Constants.hpp"
#include "hydrodynamics/Sphere.hpp"
#include "hydrodynamics/Ellipsoid.hpp"
#include "applications/hydrodynamics/CompositeShape.hpp"
#include "math/LU.hpp"
#include "utils/simError.h"



namespace OpenMD {

  /**
   * Reference:
   * Beatriz Carrasco and Jose Gracia de la Torre; "Hydrodynamic Properties of
   * Rigid Particles: Comparison of Different Modeling and Computational
   * Procedures", Biophysical Journal, 75(6), 3044, 1999   //(Hydro framework)
   *
   * Xiuquan Sun, Teng Lin, and J. Daniel Gezelter; "Langevin dynamics for rigid
   * bodies of arbitrary shape", J. Chem. Phys. 128, 234107 (2008)   //(Hydro framework)
   *
   * Beatriz Carrasco and Jose Garcia de la Torre and Peter Zipper; "Calculation
   * of hydrodynamic properties of macromolecular bead models with overlapping
   * spheres", Eur Biophys J (1999) 28: 510-515   //(overlapping beads and overlapping volume)
   *
   * http://mathworld.wolfram.com/Sphere-SphereIntersection.html   //(overlapping volume between two beads (spheres))
   *
   * Zuk, P. J., E. Wajnryb, K. A. Mizerski, and P. Szymczak; “Rotne–Prager–Yamakawa
   * Approximation for Different-Sized Particles in Application to Macromolecular
   * Bead Models.”, Journal of Fluid Mechanics, 741 (2014)   //(non-overlapping and overlapping translational-translational mobility tensors)
   *
   * Steven Harvey and Jose Garcia de la Torre; "Coordinate Systems for Modeling
   * the Hydrodynamic Resistance and Diffusion Coefficients of Irregularly Shaped
   * Rigid Macromolecules", Macromolecules 1980 13 (4), 960-964   //(center of resistance and
   * center of diffusion)
   *
   **/
  ApproximationModel::ApproximationModel(StuntDouble* sd, SimInfo* info) :
    HydrodynamicsModel(sd, info){
  }

  void ApproximationModel::init() {
    if (!createBeads(beads_)) {
      sprintf(painCave.errMsg,
              "ApproximationModel::init() : Could not create beads\n");
      painCave.isFatal = 1;
      simError();
    }
  }

  bool ApproximationModel::calcHydroProps(Shape* shape, RealType viscosity,
                                          RealType temperature) {

    HydroProp* cr = new HydroProp();
    HydroProp* cd = new HydroProp();
    HydroProp* com = new HydroProp();
    calcHydroPropsAtCRandAtCDandAtCOM(beads_, viscosity, temperature, cr, cd, com);
    setCR(cr);
    setCD(cd);
    setCOM(com);
    return true;
  }

  bool ApproximationModel::calcHydroPropsAtCRandAtCDandAtCOM(std::vector<BeadParam>& beads,
                                              RealType viscosity,
                                              RealType temperature,
                                              HydroProp* cr, HydroProp* cd, HydroProp* com) {

    std::cout << "\n";

    unsigned int nbeads = beads.size();
    DynamicRectMatrix<RealType> B(3*nbeads, 3*nbeads);
    DynamicRectMatrix<RealType> C(3*nbeads, 3*nbeads);
    Mat3x3d I;
    I(0, 0) = 1.0;
    I(1, 1) = 1.0;
    I(2, 2) = 1.0;
    RealType volume_overlap = 0.0;
    RealType overlap_beads = 0.0;
    RealType overlap_percent = 0.0;

    for (std::size_t i = 0; i < nbeads; ++i) {
      for (std::size_t j = 0; j < nbeads; ++j) {
         //checking if the beads' radii are non-negative values.
        if (beads[i].radius < 0 || beads[j].radius < 0){
          sprintf(painCave.errMsg, "There are beads with negative radius. Starting from index 0,\
 check bead (%lu) and/or bead (%lu).\n", i, j);
          painCave.isFatal = 1;
          simError();
        }
        //if the bead's radius is below 1.0e-14, substitute by 1.0e-14;
        //to avoid problem in the self-interaction part (i.e., to not divide by zero)
        if (beads[i].radius < 1.0e-14){
          beads[i].radius = 1.0e-14;
        }
        else if (beads[j].radius < 1.0e-14){
          beads[j].radius = 1.0e-14;
        }

        Mat3x3d Tij;
        if (i != j ) {   //non-self interaction: divided in overlapping and non-overlapping beads; the transitions among them are continuous
          Vector3d Rij = beads[i].pos - beads[j].pos;
          RealType rij = Rij.length();
          RealType rij2 = rij * rij;
          if (rij >= (beads[i].radius + beads[j].radius)) {      //non-overlapping beads
            RealType sumSigma2OverRij2 = ((beads[i].radius*beads[i].radius) +
                                               (beads[j].radius*beads[j].radius)) / rij2;
            Mat3x3d tmpMat;
            tmpMat = outProduct(Rij, Rij) / rij2;
            RealType constant = 8.0 * Constants::PI * viscosity * rij;
            RealType tmp1 = 1.0 + sumSigma2OverRij2/3.0;
            RealType tmp2 = 1.0 - sumSigma2OverRij2;
            Tij = (tmp1 * I + tmp2 * tmpMat ) / constant;
        }      //overlapping beads, part I
          else if ( rij > fabs(beads[i].radius - beads[j].radius) && rij < (beads[i].radius + beads[j].radius) ) {
            std::cout << "There is overlapping between beads: (" << i << ") and (" << j << ")" << std::endl;
            overlap_beads +=1;  //counting overlapping beads
            RealType sum_sigma = (beads[i].radius + beads[j].radius);
            RealType subtr_sigma = (beads[i].radius - beads[j].radius);
            RealType subtr_sigma_sqr = subtr_sigma * subtr_sigma;

            RealType constant = 6.0 * Constants::PI * viscosity * (beads[i].radius * beads[j].radius);

            RealType rij3 = rij2 * rij;

            Mat3x3d tmpMat;
            tmpMat = outProduct(Rij, Rij) / rij2;

            RealType tmp_var1 = subtr_sigma_sqr + 3.0 * rij2;
            RealType tmp1overlap = ( 16.0 * rij3 * sum_sigma - tmp_var1 * tmp_var1 )/( 32.0 * rij3 );

            RealType tmp_var2 = subtr_sigma_sqr - rij2;
            RealType tmp2overlap = ( 3.0 * tmp_var2 * tmp_var2 )/( 32.0 * rij3 );

            Tij = (tmp1overlap * I + tmp2overlap * tmpMat) / constant;

            RealType volumetmp1 = (-rij + sum_sigma) * (-rij + sum_sigma);
            RealType volumetmp2 = (rij2 + 2.0 * (rij * beads[i].radius) - 3.0 * (beads[i].radius * beads[i].radius) +
               2 * (rij * beads[j].radius) + 6 * (beads[i].radius * beads[j].radius) - 3 * (beads[j].radius * beads[j].radius));
            volume_overlap += (Constants::PI/(12.0 * rij)) * volumetmp1 * volumetmp2;
        }
          else{      //overlapping beads, part II: one bead inside the other
            if ( (beads[i].radius >= beads[j].radius) ) {
              std::cout << "Bead: (" << j << ") is inside bead (" << i << ")" << std::endl;
              overlap_beads +=1;  //counting overlapping beads
              RealType constant = 1.0 / (6.0 * Constants::PI * viscosity * beads[i].radius);
              Tij(0, 0) = constant;
              Tij(1, 1) = constant;
              Tij(2, 2) = constant;

              RealType bead_rad_cub = beads[j].radius * beads[j].radius * beads[j].radius;

              volume_overlap += (4.0/3.0) * Constants::PI * bead_rad_cub;

            } else{
              std::cout << "Bead: (" << i << ") is inside bead (" << j << ")" << std::endl;
              overlap_beads +=1;  //counting overlapping beads
              RealType constant = 1.0 / (6.0 * Constants::PI * viscosity * beads[j].radius);
              Tij(0, 0) = constant;
              Tij(1, 1) = constant;
              Tij(2, 2) = constant;

              RealType bead_rad_cub = beads[i].radius * beads[i].radius * beads[i].radius;

              volume_overlap += (4.0/3.0) * Constants::PI * bead_rad_cub;
            }
          }
        }else {   //self interaction, there is no overlapping volume (volume_overlap)
          RealType constant = 1.0 / (6.0 * Constants::PI * viscosity * beads[i].radius);
          Tij(0, 0) = constant;
          Tij(1, 1) = constant;
          Tij(2, 2) = constant;
        }
        B.setSubMatrix(i*3, j*3, Tij);
      }
    }


    //overlapping beads percentage
    //#overlap_beads counts (i,j) and (j,i) beads; and N*(N-1) counts (Tij) and (Tji) elements  (i!=j)
    if (overlap_beads > 0){
      overlap_percent = (overlap_beads * 1.0/(nbeads * (nbeads-1)))*100;  //(overlap_beads/(N*(N-1)))*100
    }

    //invert B Matrix
    invertMatrix(B, C);  //B is modified during the inversion

    //prepare U Matrix relative to arbitrary origin O(0.0, 0.0, 0.0); from the generated geometry file .xyz
    std::vector<Mat3x3d> U;
    for (unsigned int i = 0; i < nbeads; ++i) {
      Mat3x3d currU;
      currU.setupSkewMat(beads[i].pos);
      U.push_back(currU);
    }

    //calculate Xi matrix at arbitrary origin O
    Mat3x3d Xiott;
    Mat3x3d Xiorr;
    Mat3x3d Xiotr;

    //calculate the total volume, discounting the overlapping beads computed above (assuming that all overlaps are binary)
    RealType volume = -0.5 * volume_overlap;  //volume_overlap double count Vij = Vji
    for (std::vector<BeadParam>::iterator iter = beads.begin();
         iter != beads.end(); ++iter) {
      volume += 4.0/3.0 * Constants::PI * pow((*iter).radius,3);
    }

    for (std::size_t i = 0; i < nbeads; ++i) {
      for (std::size_t j = 0; j < nbeads; ++j) {
        Mat3x3d Cij;
        C.getSubMatrix(i*3, j*3, Cij);

        Xiott += Cij;
        Xiotr += U[i] * Cij;
	// Uncorrected here.  Volume correction is added after we
	// assemble Xiorr
        Xiorr += -U[i] * Cij * U[j];
      }
    }

    // Add the volume correction
    Xiorr += (RealType(6.0) * viscosity * volume) * I;

    Xiott *= Constants::viscoConvert;
    Xiotr *= Constants::viscoConvert;
    Xiorr *= Constants::viscoConvert;

    //calculate center of resistance
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

    //invert tmp Matrix
    invertMatrix(tmp, tmpInv);

    // center of resistance
    Vector3d ror = tmpInv * tmpVec;

    //calculate Resistance Tensor at center of resistance
    Mat3x3d Uor;
    Uor.setupSkewMat(ror);

    Mat3x3d Xirtt;
    Mat3x3d Xirrr;
    Mat3x3d Xirtr;

    //Resistance tensors at the center of resistance
    Xirtt = Xiott;
    Xirtr = (Xiotr - Uor * Xiott);
    Xirrr = Xiorr - Uor * Xiott * Uor + Xiotr * Uor - Uor * Xiotr.transpose();

    //calculate Diffusion tensors at center of resistance
    SquareMatrix<RealType,6> Xir6x6;
    SquareMatrix<RealType,6> Dr6x6;

    Xir6x6.setSubMatrix(0, 0, Xirtt);
    Xir6x6.setSubMatrix(0, 3, Xirtr.transpose());
    Xir6x6.setSubMatrix(3, 0, Xirtr);
    Xir6x6.setSubMatrix(3, 3, Xirrr);

    invertMatrix(Xir6x6, Dr6x6);  //Xir6x6 is modified during the inversion

    RealType kt = Constants::kb * temperature ; // kt in kcal mol^-1

    Dr6x6 *= kt;

    Mat3x3d Drtt;
    Mat3x3d Drtr;
    Mat3x3d Drrt;
    Mat3x3d Drrr;
    Dr6x6.getSubMatrix(0, 0, Drtt);
    Dr6x6.getSubMatrix(0, 3, Drrt);
    Dr6x6.getSubMatrix(3, 0, Drtr);
    Dr6x6.getSubMatrix(3, 3, Drrr);


    //Write data to .hydro file
    Mat6x6d Xi;

    cr->setCOR(ror);

    Xi.setSubMatrix(0, 0, Xirtt);
    Xi.setSubMatrix(0, 3, Xirtr.transpose());
    Xi.setSubMatrix(3, 0, Xirtr);
    Xi.setSubMatrix(3, 3, Xirrr);

    cr->setXi(Xi);

    cr->setD(Dr6x6);

    std::cout << "\n";
    std::cout << "-----------------------------------------\n";
    std::cout << "viscosity = " << viscosity << " Poise" << std::endl;
    std::cout << "temperature = " << temperature << " K" <<std::endl;
    std::cout << "-----------------------------------------\n";
    std::cout << "Percentage of overlapping beads = " << overlap_percent << " %" << std::endl;
    std::cout << "-----------------------------------------\n\n";
    std::cout << "The centers are based on the beads generated by Hydro (.xyz file), i.e.," << std::endl;
    std::cout << "not based on the geometry in the .omd file.\n" << std::endl;
    std::cout << "-----------------------------------------\n\n";
    std::cout << "Center of resistance :" << std::endl;
    std::cout << ror << "\n" << std::endl;
    std::cout << "-----------------------------------------\n\n";
    std::cout << "Resistance tensor at center of resistance\n" << std::endl;
    std::cout << "translation [kcal.fs/(mol.A^2)]:" << std::endl;
    std::cout << Xirtt << std::endl;
    std::cout << "rotation-translation [kcal.fs/(mol.A.radian)]:" << std::endl;
    std::cout << Xirtr.transpose() << std::endl;
    std::cout << "translation-rotation [kcal.fs/(mol.A.radian)]:" << std::endl;
    std::cout << Xirtr << std::endl;
    std::cout << "rotation [kcal.fs/(mol.radian^2)]:" << std::endl;
    std::cout << Xirrr << std::endl;
    std::cout << "-----------------------------------------\n\n";
    std::cout << "Diffusion tensor at center of resistance\n" << std::endl;
    std::cout << "translation (A^2 / fs):" << std::endl;
    std::cout << Drtt << std::endl;
    std::cout << "rotation-translation (A.radian / fs):" << std::endl;
    std::cout << Drrt << std::endl;
    std::cout << "translation-rotation (A.radian / fs):" << std::endl;
    std::cout << Drtr << std::endl;
    std::cout << "rotation (radian^2 / fs):" << std::endl;
    std::cout << Drrr << std::endl;
    std::cout << "-----------------------------------------\n\n";


    //calculate center of diffusion using the same arbitrary origin as above (from the generated geometry file .xyz)
    SquareMatrix<RealType,6> Xio6x6;
    SquareMatrix<RealType,6> Do6x6;

    Xio6x6.setSubMatrix(0, 0, Xiott);
    Xio6x6.setSubMatrix(0, 3, Xiotr.transpose());
    Xio6x6.setSubMatrix(3, 0, Xiotr);
    Xio6x6.setSubMatrix(3, 3, Xiorr);

    invertMatrix(Xio6x6, Do6x6);  //Xio6x6 is modified during the inversion

    // (kt) in kcal mol^-1
    Do6x6*= kt;

    Mat3x3d Dott;
    Mat3x3d Dotr;
    Mat3x3d Dort;
    Mat3x3d Dorr;
    Do6x6.getSubMatrix(0, 0, Dott);
    Do6x6.getSubMatrix(0, 3, Dort);
    Do6x6.getSubMatrix(3, 0, Dotr);
    Do6x6.getSubMatrix(3, 3, Dorr);

    //center of diffusion
    tmp(0, 0) = Dorr(1, 1) + Dorr(2, 2);
    tmp(0, 1) = -Dorr(0, 1);
    tmp(0, 2) = -Dorr(0, 2);
    tmp(1, 0) = -Dorr(0, 1);
    tmp(1, 1) = Dorr(0, 0)  + Dorr(2, 2);
    tmp(1, 2) = -Dorr(1, 2);
    tmp(2, 0) = -Dorr(0, 2);
    tmp(2, 1) = -Dorr(1, 2);
    tmp(2, 2) = Dorr(1, 1) + Dorr(0, 0);

    //Vector3d tmpVec;
    tmpVec[0] = Dotr(1, 2) - Dotr(2, 1);
    tmpVec[1] = Dotr(2, 0) - Dotr(0, 2);
    tmpVec[2] = Dotr(0, 1) - Dotr(1, 0);

    //invert tmp Matrix
    invertMatrix(tmp, tmpInv);

    //center of difussion
    Vector3d rod = tmpInv * tmpVec;

    //calculate Diffusion Tensor at center of diffusion
    Mat3x3d Uod;
    Uod.setupSkewMat(rod);

    Mat3x3d Ddtt; //translational diffusion tensor at diffusion center
    Mat3x3d Ddtr; //rotational diffusion tensor at diffusion center
    Mat3x3d Ddrr; //translation-rotation coupling diffusion tensor at diffusion tensor (which is symmetric)

    Ddtt = Dott - Uod * Dorr * Uod + Dotr.transpose() * Uod - Uod * Dotr;
    Ddrr = Dorr;
    Ddtr = Dotr + Dorr * Uod;

    SquareMatrix<RealType, 6> Dd;
    Dd.setSubMatrix(0, 0, Ddtt);
    Dd.setSubMatrix(0, 3, Ddtr.transpose());
    Dd.setSubMatrix(3, 0, Ddtr);
    Dd.setSubMatrix(3, 3, Ddrr);

    SquareMatrix<RealType, 6> Xid;
    invertMatrix(Dd, Xid);  //Dd is modified during the inversion

    // D=(kt) Z^-1  ==>  (kt) (D^-1) = Z;  t=temperature; D^-1 = Xid
    Xid *= kt;

    Mat3x3d Xidtt;
    Mat3x3d Xidrt;
    Mat3x3d Xidtr;
    Mat3x3d Xidrr;
    Xid.getSubMatrix(0, 0, Xidtt);
    Xid.getSubMatrix(0, 3, Xidrt);
    Xid.getSubMatrix(3, 0, Xidtr);
    Xid.getSubMatrix(3, 3, Xidrr);


    //Write data to .hydro file
    Mat6x6d Did;

    cd->setCOR(rod);

    cd->setXi(Xid);

    Did.setSubMatrix(0, 0, Ddtt);
    Did.setSubMatrix(0, 3, Ddtr.transpose());
    Did.setSubMatrix(3, 0, Ddtr);
    Did.setSubMatrix(3, 3, Ddrr);

    cd->setD(Did);

    std::cout << "Center of diffusion: " << std::endl;
    std::cout << rod << "\n" << std::endl;
    std::cout << "-----------------------------------------\n\n";
    std::cout << "Diffusion tensor at center of diffusion \n " << std::endl;
    std::cout << "translation (A^2 / fs) :" << std::endl;
    std::cout << Ddtt << std::endl;
    std::cout << "rotation-translation (A.radian / fs):" << std::endl;
    std::cout << Ddtr.transpose() << std::endl;
    std::cout << "translation-rotation (A.radian / fs):" << std::endl;
    std::cout << Ddtr << std::endl;
    std::cout << "rotation (radian^2 / fs):" << std::endl;
    std::cout << Ddrr << std::endl;
    std::cout << "-----------------------------------------\n\n";
    std::cout << "Resistance tensor at center of diffusion \n " << std::endl;
    std::cout << "translation [kcal.fs/(mol.A^2)]:" << std::endl;
    std::cout << Xidtt << std::endl;
    std::cout << "rotation-translation [kcal.fs/(mol.A.radian)]:" << std::endl;
    std::cout << Xidrt << std::endl;
    std::cout << "translation-rotation [kcal.fs/(mol.A.radian)]:" << std::endl;
    std::cout << Xidtr << std::endl;
    std::cout << "rotation [kcal.fs/(mol.radian^2)]:" << std::endl;
    std::cout << Xidrr << std::endl;
    std::cout << "-----------------------------------------\n\n";


    //computing center of mass using the arbitrary origin as above (from the generated .xyz geometry)
    Vector3d refMolCom;
    RealType molMass;
    molMass = 0.0;
    refMolCom=V3Zero;

    for(std::size_t i = 0; i < nbeads; ++i){
      refMolCom += beads[i].pos * beads[i].mass;
      molMass += beads[i].mass;
    }

    refMolCom /= molMass;

    Mat3x3d UoCOM;
    UoCOM.setupSkewMat(refMolCom);

    Mat3x3d Xicomtt;
    Mat3x3d Xicomrr;
    Mat3x3d Xicomtr;

    //Resistance tensors at the center of mass
    Xicomtt = Xiott;
    Xicomtr = (Xiotr - UoCOM * Xiott);
    Xicomrr = Xiorr - UoCOM * Xiott * UoCOM + Xiotr * UoCOM - UoCOM * Xiotr.transpose();

    //calculate Diffusion Tensor at center of mass
    Mat3x3d Dcomtt; //translational diffusion tensor at center of mass
    Mat3x3d Dcomtr; //rotational diffusion tensor at center of mass
    Mat3x3d Dcomrr; //translation-rotation coupling diffusion tensor at center of mass

    Dcomtt = Dott - UoCOM * Dorr * UoCOM + Dotr.transpose() * UoCOM - UoCOM * Dotr;
    Dcomrr = Dorr;
    Dcomtr = Dotr + Dorr * UoCOM;


    //Write data to .hydro file
    Mat6x6d Xcm, Dcm;

    com->setCOR(refMolCom);

    Xcm.setSubMatrix(0, 0, Xicomtt);
    Xcm.setSubMatrix(0, 3, Xicomtr.transpose());
    Xcm.setSubMatrix(3, 0, Xicomtr);
    Xcm.setSubMatrix(3, 3, Xicomrr);

    com->setXi(Xcm);

    Dcm.setSubMatrix(0, 0, Dcomtt);
    Dcm.setSubMatrix(0, 3, Dcomtr.transpose());
    Dcm.setSubMatrix(3, 0, Dcomtr);
    Dcm.setSubMatrix(3, 3, Dcomrr);

    com->setD(Dcm);


    std::cout << "Center of mass :" << std::endl;
    std::cout << refMolCom << "\n" << std::endl;
    std::cout << "-----------------------------------------\n\n";
    std::cout << "Resistance tensor at center of mass\n" << std::endl;
    std::cout << "translation [kcal.fs/(mol.A^2)]:" << std::endl;
    std::cout << Xicomtt << std::endl;
    std::cout << "rotation-translation [kcal.fs/(mol.A.radian)]:" << std::endl;
    std::cout << Xicomtr.transpose() << std::endl;
    std::cout << "translation-rotation [kcal.fs/(mol.A.radian)]:" << std::endl;
    std::cout << Xicomtr << std::endl;
    std::cout << "rotation [kcal.fs/(mol.radian^2)]:" << std::endl;
    std::cout << Xicomrr << std::endl;
    std::cout << "-----------------------------------------\n\n";
    std::cout << "Diffusion tensor at center of mass\n" << std::endl;
    std::cout << "translation (A^2 / fs):" << std::endl;
    std::cout << Dcomtt << std::endl;
    std::cout << "rotation-translation (A.radian / fs):" << std::endl;
    std::cout << Dcomtr.transpose() << std::endl;
    std::cout << "translation-rotation (A.radian / fs):" << std::endl;
    std::cout << Dcomtr << std::endl;
    std::cout << "rotation (radian^2 / fs):" << std::endl;
    std::cout << Dcomrr << std::endl;


    return true;
  }


  void ApproximationModel::writeBeads(std::ostream& os) {
    std::vector<BeadParam>::iterator iter;
    os << beads_.size() << std::endl;
    os << "Generated by Hydro" << std::endl;
    for (iter = beads_.begin(); iter != beads_.end(); ++iter) {
      os << iter->atomName << "\t"
         << iter->pos[0] << "\t"
         << iter->pos[1] << "\t"
         << iter->pos[2] << std::endl;
    }
  }
}
