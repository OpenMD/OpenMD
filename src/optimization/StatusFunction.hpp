#ifndef OPTIMIZATION_STATUSFUNCTION_HPP
#define OPTIMIZATION_STATUSFUNCTION_HPP

#include <config.h>

#include "brains/Stats.hpp"
#include "io/DumpWriter.hpp"
#include "io/StatWriter.hpp"

namespace OpenMD {
  class StatusFunction {
  public:
    virtual ~StatusFunction() {}
    virtual void writeStatus(int, int, const DynamicVector<RealType>&,
                             RealType) {}
  };

  //! No status
  class NoStatus : public StatusFunction {
  public:
    virtual void writeStatus(int, int, const DynamicVector<RealType>&,
                             RealType) {};
  };

  class DumpStatusFunction : public StatusFunction {
  public:
    DumpStatusFunction(SimInfo* info) : StatusFunction(), info_(info) {
      stats      = new Stats(info_);
      dumpWriter = new DumpWriter(info_);
      Stats::StatsBitSet mask;
      mask.set(Stats::TIME);
      mask.set(Stats::POTENTIAL_ENERGY);
      stats->setStatsMask(mask);
      statWriter = new StatWriter(info_->getStatFileName(), stats);
    }

    virtual void writeStatus(int functionCount, int,
                             const DynamicVector<RealType>&, RealType) {
      Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
      curSnapshot->setTime(functionCount);
      stats->collectStats();
      statWriter->writeStat();

      dumpWriter->writeDumpAndEor();
    }
    ~DumpStatusFunction() {
      delete stats;
      delete dumpWriter;
      delete statWriter;
    }

  private:
    SimInfo* info_ {nullptr};
    Stats* stats;
    DumpWriter* dumpWriter;
    StatWriter* statWriter;
  };
}  // namespace OpenMD

#endif
