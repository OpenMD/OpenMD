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
#ifndef APPLICATIONS_DYNAMICPROPS_CORRELATIONFUNCTION_HPP
#define APPLICATIONS_DYNAMICPROPS_CORRELATIONFUNCTION_HPP

#include <string>
#include <vector>

#include "brains/SimInfo.hpp"
#include "brains/BlockSnapshotManager.hpp"

#include "primitives/StuntDouble.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"

namespace oopse {

/**
 * @class CorrelationFunction CorrelationFunction.hpp "applications/dynamicProps/CorrelationFunction"
 * @brief Base class for Correlation function
 */
 
class CorrelationFunction {
    public:
        CorrelationFunction(SimInfo* info, const std::string& filename, 
            const std::string& sele1, const std::string& sele2, int storageLayout);
        
        void doCorrelate();


        void setOutputName(const std::string& filename) {
            outputFilename_ = filename;
        }

        const std::string& getOutputFileName() const {
            return outputFilename_;
        }


        const std::string& getCorrFuncType() const {
            return corrFuncType_;
        }

        void setCorrFuncType(const std::string& type) {
            corrFuncType_ = type;
        }

        void setExtraInfo(const std::string& extra) {
            extra_ = extra;
        }
            
    protected:
        
        virtual void preCorrelate();        
        virtual void postCorrelate();

        double deltaTime_;
        int nTimeBins_;
        std::vector<double> histogram_;
        std::vector<int> count_;
        std::vector<int> time_;
        
        SimInfo* info_;
        int storageLayout_;
        std::string dumpFilename_;        
        SelectionManager seleMan1_;
        SelectionManager seleMan2_;          
        
    private:

        void correlateBlocks(int block1, int block2);
        void correlateFrames(int frame1, int frame2);
        void updateFrame(int frame);
        
        virtual void writeCorrelate();

        virtual double calcCorrVal(StuntDouble* sd1, int frame1, StuntDouble* sd2, int frame2) = 0;

        virtual void validateSelection(const SelectionManager& seleMan) {}


        std::string selectionScript1_;
        std::string selectionScript2_;
        
        SelectionEvaluator evaluator1_;
        SelectionEvaluator evaluator2_;

        BlockSnapshotManager* bsMan_;        

        std::string outputFilename_;

        std::string corrFuncType_;
        std::string extra_;
};

}
#endif
