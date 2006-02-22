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

#include "applications/hydrodynamics/MoleculeShape.hpp" 
#include "utils/MemoryUtils.hpp"

namespace oopse {

Spheric::Spheric(Vector3d origin, double radius) : origin_(origin), radius_(radius){

}

bool Spheric::isInterior(Vector3d pos) {
    Vector3d r = pos - origin_;

    bool result;
    if (r.length() < radius_)
        result = true;
    else
        result = false;
    
    return result;
}

std::pair<Vector3d, Vector3d> Spheric::getBox() {
    std::pair<Vector3d, Vector3d>  boundary;
    Vector3d r(radius_, radius_, radius_);
    boundary.first = origin_ - r;
    boundary.first = origin_ + r;
    return boundary;
}
Ellipsoid::Ellipsoid(Vector3d origin, double radius, double ratio, Mat3x3d rotMat) 
    : origin_(origin), a_(radius), b_(radius*ratio), rotMat_(rotMat) {

}
bool Ellipsoid::isInterior(Vector3d pos) {
    Vector3d r = pos - origin_;
    Vector3d rbody = rotMat_ * r;
    double xovera = rbody[0]/a_;
    double yovera = rbody[1]/a_;
    double zoverb = rbody[2]/b_;

    bool result;
    if (xovera*xovera + yovera*yovera + zoverb*zoverb < 1)
        result = true;
    else
        result = false;

    return result;
        
}

std::pair<Vector3d, Vector3d> Ellipsoid::getBox() {

    std::pair<Vector3d, Vector3d>  boundary;
    //make a cubic box
    double rad  = a_ > b_ ? a_ : b_; 
    Vector3d r(rad, rad, rad);
    boundary.first = origin_ - r;
    boundary.first = origin_ + r;
    return boundary;
}


MoleculeShape::MoleculeShape(Molecule* mol) {
    std::vector<Atom*>::iterator ai; 
    Atom* atom;
    for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
        AtomType* atomType = atom->getAtomType();
        Shape* currShape = NULL;
        if (atomType->isGayBerne()) {
            DirectionalAtomType* dAtomType = dynamic_cast<DirectionalAtomType*>(atomType);

            GenericData* data = dAtomType->getPropertyByName("GayBerne");
            if (data != NULL) {
                GayBerneParamGenericData* gayBerneData = dynamic_cast<GayBerneParamGenericData*>(data);

                if (gayBerneData != NULL) {  
                    GayBerneParam gayBerneParam = gayBerneData->getData();
                    currShape = new Ellipsoid(atom->getPos(), gayBerneParam.GB_sigma, gayBerneParam.GB_l2b_ratio, atom->getA());
                } else {
                    sprintf( painCave.errMsg,
                           "Can not cast GenericData to GayBerneParam\n");
                    painCave.severity = OOPSE_ERROR;
                    painCave.isFatal = 1;
                    simError();   
                }
            } else {
        	  sprintf( painCave.errMsg, "Can not find Parameters for GayBerne\n");
        	  painCave.severity = OOPSE_ERROR;
        	  painCave.isFatal = 1;
        	  simError();    
            }            
        } else if (atomType->isLennardJones()){
            GenericData* data = atomType->getPropertyByName("LennardJones");
            if (data != NULL) {
                LJParamGenericData* ljData = dynamic_cast<LJParamGenericData*>(data);

                if (ljData != NULL) {
                    LJParam ljParam = ljData->getData();
                    currShape = new Spheric(atom->getPos(), ljParam.sigma/2.0);
            } else {
                sprintf( painCave.errMsg,
                "Can not cast GenericData to LJParam\n");
                painCave.severity = OOPSE_ERROR;
                painCave.isFatal = 1;
                simError();          
                }       
            }

        }

        if (currShape != NULL)
            shapes_.push_back(currShape);

    }

}

MoleculeShape::~MoleculeShape() {
    MemoryUtils::deletePointers(shapes_);
}
bool MoleculeShape::isInterior(Vector3d pos) {
    bool result = false;
    std::vector<Shape*>::iterator iter;
    for (iter = shapes_.begin(); iter != shapes_.end(); ++ iter) {
        if ((*iter)->isInterior(pos)) {
            result = true;
            break;
        }
    }

    return result;
}

template<class Cont, class Predict>
void swap_if(Cont& b1, Cont& b2, Predict predict) {
    unsigned int size = b1.size();
    assert(size == b2.size());
    for (unsigned int i = 0 ; i < size; ++i) {
        if (predict(b1[i], b2[i]))
            std::swap(b1[i], b2[i]);
    }

}

std::pair<Vector3d, Vector3d> MoleculeShape::getBox() {
    std::vector<Shape*>::iterator iter = shapes_.begin();
    std::pair<Vector3d, Vector3d>  boundary = (*iter)->getBox();
    for (++iter; iter != shapes_.end(); ++iter) {
        std::pair<Vector3d, Vector3d> currBoundary = (*iter)->getBox();
            swap_if(boundary.first, currBoundary.first, std::less<double>());
            swap_if(boundary.second, currBoundary.second, std::greater<double>());        
    }

    return boundary;
}


}
