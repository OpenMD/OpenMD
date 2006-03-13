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

#ifndef APPLICATION_HYDRODYNAMICS_MOLECULESHAPE_HPP
#define APPLICATION_HYDRODYNAMICS_MOLECULESHAPE_HPP

#include "math/Vector3.hpp"
#include "primitives/Molecule.hpp"
namespace oopse {

/** @class Shape*/
class Shape {
    public:
        virtual bool isInterior(Vector3d pos) = 0;
        virtual std::pair<Vector3d, Vector3d> getBox() = 0;
        virtual bool hasAnalyticalSolution() = 0;
        //virtual Mat3x3d getTranslationalFriction() = 0;
        //virtual Mat3x3d getRotationalFriction() = 0 ;
};

/** @class Spheric */
class Spheric : public Shape {
    public:
        Spheric(Vector3d origin, double radius);
        virtual bool isInterior(Vector3d pos);
        virtual std::pair<Vector3d, Vector3d> getBox();
        virtual bool hasAnalyticalSolution() {return true;}
        //virtual Mat3x3d getTranslationalFriction();
        //virtual Mat3x3d getRotationalFriction();

    private:
        Vector3d origin_;
        double radius_;
};

/** @class Ellipsoid */
class Ellipsoid : public Shape{
    public:
        Ellipsoid(Vector3d origin, double radius, double ratio, Mat3x3d rotMat);
        virtual bool isInterior(Vector3d pos);
        virtual std::pair<Vector3d, Vector3d> getBox();
        virtual bool hasAnalyticalSolution() {return true;}
        //virtual Mat3x3d getTranslationalFriction() {
        //    Mat3x3d Xirr;
        //}
        
        //virtual Mat3x3d getFrictionCoff() {

        //}

     private:


        /**
         * calculate the ratio of friction coeffiction constant between ellipsoid and spheric 
         * with same volume.
         * @param m
         * @param n 
         * @note 
         * Reference:
         *
         * (1) Victor A. Bloomfield, On-Line Biophysics Textbook, Volume: Separations and Hydrodynamics
         * Chapter 1,Survey of Biomolecular Hydrodynamics
         * http://www.biophysics.org/education/vbloomfield.pdf 
         * (2) F. Perrin , J. Phys. Radium, [7] 5, 497-511, 1934
         * (3) F. Perrin, J. Phys. Radium, [7] 7, 1-11, 1936
         */
        void calcFrictionRatio(double m, double n, double& Ft, double& Frm, double& Frn) {
            double q = n/m;
            if (q > 1.0) {//prolate
                Ft = sqrt(1-q*q)/(pow(q, 2.0/3.0)*log((1 + sqrt(1-q*q))/q));
                Frm = 4*(1-q*q)/(3*(2 - 2*pow(q, 4.0/3.0)/Ft)); //not sure
                Frn = 4*(1-q*q*q*q) /(3*q*q*(2*pow(q, -2.0/3.0)*(2-q*q)/Ft-2));
            } else {//oblate
                Ft = sqrt(1-q*q)/(pow(q, 2.0/3.0)*atan(sqrt(q*q-1)));
                Frm = 4*(1-q*q)/(3*(2 - 2*pow(q, 4.0/3.0)/Ft)); //not sure
                Frn = 4*(1-q*q*q*q) /(3*q*q*(2*pow(q, -2.0/3.0)*(2-q*q)/Ft-2));
            }
        }
        Vector3d origin_;
        double a_;
        double b_;
        Mat3x3d rotMat_;
};

/** @class StuntDoubleShape */
class StuntDoubleShape {
    public:
        StuntDoubleShape(StuntDouble* sd);
        ~StuntDoubleShape();
        bool hasAnalyticalSolution();
        std::pair<Vector3d, Vector3d> getBox();
        bool isInterior(Vector3d pos);
    private:
        Shape* createShape(Atom* atom);
        std::vector<Shape*> shapes_;
};


class ShapeBuilder {
    public:
        
};


}

#endif
