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

#include "applications/dynamicProps/CorrelationFunction.hpp"
#include "utils/simError.h"

namespace oopse {

CorrelationFunction::CorrelationFunction(SimInfo* info, const std::string& filename, 
    const std::string& sele1, const std::string& sele2, int storageLayout)
    : info_(info), storageLayout_(storageLayout), dumpFilename_(filename), selectionScript1_(sele1),
      selectionScript2_(sele2), evaluator1_(info), evaluator2_(info), seleMan1_(info), seleMan2_(info){

    int nAtoms = info->getNGlobalAtoms();
    int nRigidBodies = info->getNGlobalRigidBodies();
    int nStuntDoubles =   nAtoms + nRigidBodies;

    std::set<AtomType*> atomTypes = info->getUniqueAtomTypes();
    std::set<AtomType*>::iterator i;
    bool hasDirectionalAtom = false;
    bool hasMultipole = false;    
    for (i = atomTypes.begin(); i != atomTypes.end(); ++i) {
        if ((*i)->isDirectional()){
            hasDirectionalAtom = true;
        }
        if ((*i)->isMultipole()){
            hasDirectionalAtom = true;
        }
    }
    
    if (nRigidBodies > 0 || hasDirectionalAtom) {
        storageLayout_ |= DataStorage::dslAmat;
    }
    if (hasMultipole) {
        storageLayout_ |= DataStorage::dslElectroFrame;
    }
        
    bsMan_ = new BlockSnapshotMananger(info, dumpFilename_, storageLayout_);
    info->setSnapshotManager(bsMan_);

    //if selection is static, we only need to evaluate it once
    if (!evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
        validateSelection(seleMan1_);
    }
    if (!evaluator2_.isDynamic()) {
        seleMan2_.setSelectionSet(evaluator2_.evaluate());
        validateSelection(seleMan2_);
    }    

    int nframes =  bsMan_->getNFrames();
    nTimeBins_ = nframes;  /**@todo Fixed Me */
    histogram_.resize(nTimeBins_);
    count_.resize(nTimeBins_);

    //deltaTime_ = info_->;
}


void CorrelationFunction::doCorrelate() {
    preCorrelate();

    int nblocks = bsMan_->getNBlocks();

    for (int i = 0; i < nblocks; ++i) {
        bsMan_->loadBlock(i);
        for (int j = i; j < nblocks; ++j) {
            bsMan_->loadBlock(j);
            correlateBlocks(i, j);

            if (i != j) {
                //if i equals to j, the block is still used, should not unload it
                bsMan_->unloadBlock(j);
            }
        }
        bsMan_->unloadBlock(i);
    }
    
    postCorrelate();

    writeCorrelate();
}

void CorrelationFunction::correlateBlocks(int block1, int block2) {

    assert(bsMan_->isBlockActive(block1) && bsMan_->isBlockActive(block2));

    SnapshotBlock snapshotBlock1 = bsMan_->getSnapshotBlock(block1);
    SnapshotBlock snapshotBlock2 = bsMan_->getSnapshotBlock(block2);
    
    for (int i  = snapshotBlock1.first; i < snapshotBlock1.second; ++i) {

        //evaluate selection if it is dynamic
        if (evaluator1_.isDynamic()) {
            seleMan1_.setSelectionSet(evaluator1_.evaluate());
            validateSelection(seleMan1_);
        }

        //update the position or velocity of the atoms belong to rigid bodies
        updateFrame(i);
        for(int j  = snapshotBlock1.first; j < snapshotBlock1.second; ++j) {
            //evaluate selection
            if (!evaluator2_.isDynamic()) {
                seleMan2_.setSelectionSet(evaluator2_.evaluate());
                validateSelection(seleMan2_);
            }   
            //update the position or velocity of the atoms belong to rigid bodies
            updateFrame(j);

            correlateFrames(i, j);
        }
    }
}


void CorrelationFunction::correlateFrames(int frame1, int frame2) {
    int count = seleMan1_.getSelectionCount();
    assert(  count == seleMan2_.getSelectionCount());

    Snapshot* snapshot1 = bsMan_->getSnapshot(frame1);
    Snapshot* snapshot2 = bsMan_->getSnapshot(frame2);
    assert(snapshot1 && snapshot2);

    double time1 = snapshot1->getTime();
    double time2 = snapshot2->getTime();
    int timeBin = int ((time2 - time1) /deltaTime_ + 0.5);
    count_[timeBin] += count;    
    
    int i;
    int j;
    StuntDouble* sd1;
    StuntDouble* sd2;
    for (sd1 = seleMan1_.beginSelected(i), (sd2 = seleMan2_.beginSelected(j); sd1 != NULL && sd2 != NULL;
        sd1 = seleMan1_.nextSelected(i) sd2 = seleMan2_.nextSelected(j)) {
        double corrVal = calcCorrVal(sd1, frame1, sd2, frame2);
        histogram_[timeBin] += corrVal;    
    }

    //for (sd1 = seleMan1_.beginSelected(i); sd1 != NULL; sd1 = seleMan1_.nextSelected(i)) {
    //
    //    for (sd2 = seleMan2_.beginSelected(j); sd2 != NULL; sd2 = seleMan2_.nextSelected(j)) {
    //        double corrVal = calcCorrVal(sd1, frame1, sd2, frame2);
    //    }            
    //}

}


void CorrelationFunction::updateFrame(int frame){
    Molecule* mol;
    RigidBody* rb;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;

    /** @todo need improvement */    
    if (storageLayout_ & DataStorage::dslPosition) {
        for (mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)) {

            //change the positions of atoms which belong to the rigidbodies
            for (rb = mol->beginRigidBody(rbIter); rb != NULL; rb = mol->nextRigidBody(rbIter)) {
                rb->updateAtoms(frame);
            }
        }        
    }

    if (storageLayout_ & DataStorage::dslVelocity) {
        for (mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)) {

            //change the positions of atoms which belong to the rigidbodies
            for (rb = mol->beginRigidBody(rbIter); rb != NULL; rb = mol->nextRigidBody(rbIter)) {
                rb->updateAtomVel(frame);
            }
        }
        
    }
    
}


void CorrelationFunction::preCorrelate() {
    std::fill(histogram_.begin(), histogram_.end(), 0.0);
    std::fill(count_.begin(), count_.end(), 0);
}

void CorrelationFunction::postCorrelate() {
    for (int i =0 ; i < nTimeBins_; ++i) {
        if (count_[i] > 0) {
            histogram_[i] / = count_[i];
        }
    }
}


void CorrelationFunction::writeCorrelate() {
    std::ofstream ofs(outputFilename_.c_str());

    if (ofs.is_open()) {

        ofs << "#" << getCorrFuncType() << "\n";
        ofs << "#selection script1: \"" << selectionScript1_ <<"\"\tselection script2: \"" selectionScript2_ << "\"\n";
        ofs << "#extra information: " << extra_ << "\n";
        ofs << "#time\tcorrVal\n";

        for (int i = 0; i < nTimeBins_; ++i) {
            ofs << time_[i] << "\t" << histogram_[i] << "\n";
        }
            
    } else {
            sprintf(painCave.errMsg,
                "CorrelationFunction::writeCorrelate Error: fail to open %s\n", outputFilename_.c_str());
            painCave.isFatal = 1;
            simError();        
    }

    ofs.close();    
}

}
