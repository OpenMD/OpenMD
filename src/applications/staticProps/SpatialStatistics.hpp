/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [4] , Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). *
 *  Created by Charles F. Vardeman II on 11/26/05.
 *  @author  Charles F. Vardeman II 
 *  @version $Id: SlabStatistics.hpp 1850 2013-02-20 15:39:39Z gezelter $
 *
 */
#ifndef APPLICATIONS_STATICPROPS_SPATIALSTATISTICS_HPP
#define APPLICATIONS_STATICPROPS_SPATIALSTATISTICS_HPP

#include <string>
#include <vector>
#include <fstream>

#include "applications/staticProps/StaticAnalyser.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "utils/Accumulator.hpp"

using namespace std;
namespace OpenMD {
  
  enum OutputDataType {
    odtReal,
    odtVector3,
    odtArray2d,
    odtUnknownDataType
  };

  enum OutputDataHandling {
    odhAverage,
    odhTotal,
    odhLastValue,
    odhUnknownDataHandling
  };

  struct OutputData {
    string title;
    string units;
    OutputDataType dataType;
    OutputDataHandling dataHandling;
    vector<BaseAccumulator*> accumulator;
    vector<vector<BaseAccumulator*> > accumulatorArray2d;
  };

  class SpatialStatistics : public StaticAnalyser {
    
  public:
    SpatialStatistics(SimInfo* info, const string& filename,
                      const string& sele, int nbins);    
    ~SpatialStatistics();

    void addOutputData(OutputData* dat) {data_.push_back(dat);}
    virtual void process();
    virtual void processFrame(int frame);
    virtual int getBin(Vector3d pos)=0;
    virtual void processStuntDouble(StuntDouble* sd, int bin)=0;
    virtual void writeOutput();
    
  protected:
    OutputData* beginOutputData(vector<OutputData*>::iterator& i);
    OutputData* nextOutputData(vector<OutputData*>::iterator& i);
    void writeData(ostream& os, OutputData* dat, unsigned int bin);
    void writeErrorBars(ostream& os, OutputData* dat, unsigned int bin);

    Snapshot* currentSnapshot_;    
    int nProcessed_;
    string selectionScript_;
    SelectionEvaluator evaluator_;
    SelectionManager seleMan_;    
    int nBins_; 
    OutputData* counts_;
    vector<OutputData*> data_;
  };
  
  class SlabStatistics : public SpatialStatistics {
  public: 
    SlabStatistics(SimInfo* info, const string& filename,
                   const string& sele, int nbins);
    virtual ~SlabStatistics();

    virtual int getBin(Vector3d pos);
    virtual void processFrame(int frame);
  protected:
    OutputData* z_;
    Mat3x3d hmat_;
    RealType volume_;

  };

  class ShellStatistics : public SpatialStatistics {
  public: 
    ShellStatistics(SimInfo* info, const string& filename, const string& sele,
                    int nbins);
    virtual ~ShellStatistics();
    virtual int getBin(Vector3d pos);

    void setCoordinateOrigin(Vector3d co) { coordinateOrigin_ = co; }
    void setBinWidth(RealType bw) { binWidth_ = bw; }

  protected:
    OutputData* r_;
    Vector3d coordinateOrigin_;
    RealType binWidth_;
  };
  
}
#endif



