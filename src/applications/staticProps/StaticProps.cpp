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
 
#include <iostream>
#include <fstream>
#include <string>

#include "brains/Register.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "io/DumpReader.hpp"
#include "utils/simError.h"

#include "applications/staticProps/StaticPropsCmd.h"
#include "applications/staticProps/GofR.hpp"
#include "applications/staticProps/GofRAngle.hpp"
#include "applications/staticProps/GofAngle2.hpp"
#include "applications/staticProps/GofXyz.hpp"

using namespace oopse;

int main(int argc, char* argv[]){
  
    //register force fields
    registerForceFields();

    gengetopt_args_info args_info;

    //parse the command line option
    if (cmdline_parser (argc, argv, &args_info) != 0) {
        exit(1) ;
    }


    //get the dumpfile name and meta-data file name
    std::string dumpFileName = args_info.input_arg;

    std::string mdFileName = dumpFileName.substr(0, dumpFileName.rfind(".")) + ".md";

    
    std::string sele1;
    std::string sele2;

    if (args_info.sele1_given) {
        sele1 = args_info.sele1_arg;
    }else {
        char*  sele1Env= getenv("OOPSE_SELE1");
        if (sele1Env) {
            sele1 = sele1Env;
        }else {
            sprintf( painCave.errMsg,
               "neither --sele1 option nor $OOPSE_SELE1 is set");
            painCave.severity = OOPSE_ERROR;
            painCave.isFatal = 1;
            simError();
        }
    }
    
    if (args_info.sele2_given) {
        sele2 = args_info.sele2_arg;
    }else {
        char* sele2Env = getenv("OOPSE_SELE2");
        if (sele2Env) {
            sele2 = sele2Env;            
        } else {
            sprintf( painCave.errMsg,
               "neither --sele2 option nor $OOPSE_SELE2 is set");
            painCave.severity = OOPSE_ERROR;
            painCave.isFatal = 1;
            simError();        
        }
    }

    //parse md file and set up the system
    SimCreator creator;
    SimInfo* info = creator.createSim(mdFileName);

    double maxLen;
    if (args_info.length_given) {
        maxLen = args_info.length_arg;
    } else {
        Mat3x3d hmat = info->getSnapshotManager()->getCurrentSnapshot()->getHmat();
        maxLen = std::min(std::min(hmat(0, 0), hmat(1, 1)), hmat(2, 2)) /2 ;
    }
    

    RadialDistrFunc* rdf;
    if (args_info.gofr_given){
        GofR* r = new GofR(info, dumpFileName, sele1, sele2);
        
        r->setNRBins(args_info.nrbins_arg);            
        r->setLength(maxLen);
        
        rdf = r;
    } else if (args_info.r_theta_given) {
        GofRTheta* rTheta = new GofRTheta(info, dumpFileName, sele1, sele2);
          
        rTheta->setNRBins(args_info.nrbins_arg);            
        rTheta->setLength(maxLen);       
        rTheta->setNAngleBins(args_info.nanglebins_arg);

        
        rdf = rTheta;
    }
    else if (args_info.r_omega_given) {
        GofROmega* rOmega = new GofROmega(info, dumpFileName, sele1, sele2);

       
        rOmega->setNRBins(args_info.nrbins_arg);            
        rOmega->setLength(maxLen);
        rOmega->setNAngleBins(args_info.nanglebins_arg);

        rdf = rOmega;    
    } else if (args_info.theta_omega_given) {
        GofAngle2* rAngle2 = new GofAngle2(info, dumpFileName, sele1, sele2);
        rAngle2->setNAngleBins(args_info.nanglebins_arg);

        rdf = rAngle2;  
    } else if (args_info.xyz_given) {

        GofXyz* xyz = new GofXyz(info, dumpFileName, sele1, sele2);
          
        xyz->setNRBins(args_info.nrbins_arg);            
        xyz->setLength(maxLen);

        
        rdf = xyz;
    }
    

    if (args_info.output_given) {
        rdf->setOutputName(args_info.output_arg);
    }

    if (args_info.step_given) {
        rdf->setStep(args_info.step_arg);
    }

    rdf->process();

    delete rdf;    
    delete info;

    return 0;   
}

