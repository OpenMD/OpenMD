#ifndef OPTIMIZATION_STATUSFUNCTION_HPP
#define OPTIMIZATION_STATUSFUNCTION_HPP
#include "config.h"
#include "io/DumpWriter.hpp"
#include "io/StatWriter.hpp"

namespace OpenMD {
  class StatusFunction {
  public:
    virtual ~StatusFunction() {}
    virtual void writeStatus(const DynamicVector<RealType>& currentValue) { std::cerr << "doing status\n"; }    
  };

  //! No status
  class NoStatus : public StatusFunction {
  public:
    virtual void writeStatus(const DynamicVector<RealType>& currentValue) {};
  };

  class DumpStatusFunction : public StatusFunction {

  public:
    DumpStatusFunction(SimInfo* info) : StatusFunction(), info_(info) {
      dumpWriter = new DumpWriter(info_);     
      StatsBitSet mask;
      mask.set(Stats::TIME);
      mask.set(Stats::POTENTIAL_ENERGY);
      statWriter = new StatWriter(info_->getStatFileName(), mask);
    }
    virtual void writeStatus(const DynamicVector<RealType>& currentValue) {
      Snapshot* curSnapshot =info_->getSnapshotManager()->getCurrentSnapshot();
      dumpWriter->writeDumpAndEor();
      statWriter->writeStat(curSnapshot->statData);
    }
    
  private:
    SimInfo* info_;
    DumpWriter* dumpWriter;
    StatWriter* statWriter;
  };
                             
  
}
#endif
