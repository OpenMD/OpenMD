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

#include <algorithm>
#include <fstream>
#include "applications/staticProps/GofXyz.hpp"
#include "utils/simError.h"
#include "primitives/Molecule.hpp"
namespace oopse {

GofXyz::GofXyz(SimInfo* info, const std::string& filename, const std::string& sele1, const std::string& sele2, double len, int nrbins)
    : RadialDistrFunc(info, filename, sele1, sele2), len_(len), nRBins_(nrbins) {
    setOutputName(getPrefix(filename) + ".gxyz");

    deltaR_ = len_ / nRBins_;
    
    histogram_.resize(nRBins_);
    for (int i = 0 ; i < nRBins_; ++i) {
        histogram_[i].resize(nRBins_);
        for(int j = 0; j < nRBins_; ++j) {
            histogram_[i][j].resize(nRBins_);
        }
    }   

    //create atom2Mol mapping (should be other class' responsibility)
    atom2Mol_.insert(atom2Mol_.begin(), info_->getNGlobalAtoms() + info_->getNGlobalRigidBodies(), static_cast<Molecule*>(NULL));
    
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::AtomIterator ai;
    Atom* atom;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;
    
    for (mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)) {
        
        for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
            atom2Mol_[atom->getGlobalIndex()] = mol;
        }

        for (rb = mol->beginRigidBody(rbIter); rb != NULL; rb = mol->nextRigidBody(rbIter)) {
            atom2Mol_[rb->getGlobalIndex()] = mol;
        }
        
    }       
}


void GofXyz::preProcess() {
    for (int i = 0 ; i < nRBins_; ++i) {
        histogram_[i].resize(nRBins_);
        for(int j = 0; j < nRBins_; ++j) {
            std::fill(histogram_[i][j].begin(), histogram_[i][j].end(), 0);
        }
    }   
}


void GofXyz::initalizeHistogram() {
    //calculate the center of mass of the molecule of selected stuntdouble in selection1

    //determine the new coordinate set of selection1
    //v1 = Rs1 -Rcom, 
    //z = Rs1.dipole
    //x = v1 X z
    //y = z X x 
    coorSets_.clear();

    int i;
    StuntDouble* sd;
    for (sd = seleMan1_.beginSelected(i); sd != NULL; sd = seleMan1_.nextSelected(i)) {
        Vector3d rcom = getMolCom(sd);
        Vector3d rs1 = sd->getPos();
        Vector3d v1 =  rcom - rs1;
        CoorSet currCoorSet;
        currCoorSet.zaxis = sd->getElectroFrame().getColumn(2);
        v1.normalize();
        currCoorSet.zaxis.normalize();
        currCoorSet.xaxis = cross(v1, currCoorSet.zaxis);
        currCoorSet.yaxis = cross(currCoorSet.zaxis, currCoorSet.xaxis);
        coorSets_.insert(std::map<int, CoorSet>::value_type(sd->getGlobalIndex(), currCoorSet));
    }

}

void GofXyz::collectHistogram(StuntDouble* sd1, StuntDouble* sd2) {

    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12 = pos2 - pos1;
    currentSnapshot_->wrapVector(r12);

    std::map<int, CoorSet>::iterator i = coorSets_.find(sd1->getGlobalIndex());
    assert(i != coorSets_.end());
    
    double x = dot(r12, i->second.xaxis);
    double y = dot(r12, i->second.yaxis);
    double z = dot(r12, i->second.zaxis);

    int xbin = x / deltaR_;
    int ybin = y / deltaR_;
    int zbin = z / deltaR_;

    if (xbin < nRBins_ && ybin < nRBins_ && zbin < nRBins_) {
        ++histogram_[x][y][z];
    }
    
}

void GofXyz::writeRdf() {
    std::ofstream rdfStream(outputFilename_.c_str(), std::ios::binary);
    if (rdfStream.is_open()) {
        //rdfStream << "#g(x, y, z)\n";
        //rdfStream << "#selection1: (" << selectionScript1_ << ")\t";
        //rdfStream << "selection2: (" << selectionScript2_ << ")\n";
        //rdfStream << "#nRBins = " << nRBins_ << "\t maxLen = " << len_ << "deltaR = " << deltaR_ <<"\n";
        for (int i = 0; i < histogram_.size(); ++i) {
 
            for(int j = 0; j < histogram_[i].size(); ++j) {
 
                for(int k = 0;k < histogram_[i].size(); ++k) {
                    rdfStream.write(reinterpret_cast<char *>(&histogram_[i][j][k] ), sizeof(histogram_[i][j][k] ));
                }
            }
        }
        
    } else {

        sprintf(painCave.errMsg, "GofXyz: unable to open %s\n", outputFilename_.c_str());
        painCave.isFatal = 1;
        simError();  
    }

    rdfStream.close();
}

Vector3d GofXyz::getMolCom(StuntDouble* sd){
    Molecule* mol = atom2Mol_[sd->getGlobalIndex()];
    assert(mol);
    return mol->getCom();
}

}
