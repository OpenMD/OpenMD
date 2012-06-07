#ifndef OPTIMIZATION_STATUSFUNCTION_HPP
#define OPTIMIZATION_STATUSFUNCTION_HPP
#include "config.h"
#include "io/DumpWriter.hpp"
#include "io/StatWriter.hpp"

namespace OpenMD {
  class StatusFunction {
  public:
    virtual ~StatusFunction() {}
    virtual void writeStatus(int functionCount, int gradientCount, const DynamicVector<RealType>& x, RealType f) { std::cerr << "doing status\n"; }    
  };

  //! No status
  class NoStatus : public StatusFunction {
  public:
    virtual void writeStatus(int functionCount, int gradientCount, const DynamicVector<RealType>& x, RealType f) {};
  };

  class DumpStatusFunction : public StatusFunction {

  public:
    DumpStatusFunction(SimInfo* info) : StatusFunction(), info_(info), thermo(info) {
      dumpWriter = new DumpWriter(info_);     
      StatsBitSet mask;
      mask.set(Stats::TIME);
      mask.set(Stats::POTENTIAL_ENERGY);
      statWriter = new StatWriter(info_->getStatFileName(), mask);
    }
    virtual void writeStatus(int functionCount, int gradientCount, const DynamicVector<RealType>& x, RealType f) {
      Snapshot* curSnapshot =info_->getSnapshotManager()->getCurrentSnapshot();
      thermo.saveStat();
      curSnapshot->setTime(functionCount);         
      dumpWriter->writeDumpAndEor();
      statWriter->writeStat(curSnapshot->statData);
    }
    ~DumpStatusFunction() {
      delete dumpWriter;
      delete statWriter;
    }
    
  private:
    SimInfo* info_;
    DumpWriter* dumpWriter;
    StatWriter* statWriter;
    Thermo thermo;
  };
                             
  
}
#endif
