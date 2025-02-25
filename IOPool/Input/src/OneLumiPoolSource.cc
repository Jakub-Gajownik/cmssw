/*----------------------------------------------------------------------
----------------------------------------------------------------------*/
#include "PoolSource.h"
#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "FWCore/Framework/interface/LuminosityBlockPrincipal.h"
#include "FWCore/Framework/interface/EventPrincipal.h"
#include "FWCore/Framework/interface/ProductResolversFactory.h"

#include "DataFormats/Provenance/interface/ProductRegistry.h"
#include "FWCore/Framework/interface/HistoryAppender.h"

namespace edm {

  class OneLumiPoolSource : public PoolSource {
  public:
    explicit OneLumiPoolSource(ParameterSet const& pset, InputSourceDescription const& desc);

  private:
    ItemTypeInfo getNextItemType() override;
    std::shared_ptr<LuminosityBlockAuxiliary> readLuminosityBlockAuxiliary_() override;

    void readEvent_(EventPrincipal& eventPrincipal) override {
      PoolSource::readEvent_(eventPrincipal);
      eventPrincipal.setRunAndLumiNumber(eventPrincipal.run(), 1);
    }

    bool seenFirstLumi_ = false;
  };

  OneLumiPoolSource::OneLumiPoolSource(ParameterSet const& pset, InputSourceDescription const& desc)
      : PoolSource(pset, desc) {}

  std::shared_ptr<LuminosityBlockAuxiliary> OneLumiPoolSource::readLuminosityBlockAuxiliary_() {
    auto ret = PoolSource::readLuminosityBlockAuxiliary_();
    auto hist = ret->processHistoryID();
    *ret = LuminosityBlockAuxiliary(ret->run(), 1, ret->beginTime(), ret->endTime());
    ret->setProcessHistoryID(hist);
    return ret;
  }

  InputSource::ItemTypeInfo OneLumiPoolSource::getNextItemType() {
    auto type = PoolSource::getNextItemType();
    if (type == ItemType::IsLumi) {
      if (seenFirstLumi_) {
        do {
          edm::HistoryAppender historyAppender;
          auto prodReg = std::make_shared<edm::ProductRegistry>();
          prodReg->setFrozen();
          edm::ProcessConfiguration procConfig;

          LuminosityBlockPrincipal temp(
              prodReg, edm::productResolversFactory::makePrimary, procConfig, &historyAppender, 0);
          readLuminosityBlock_(temp);
          type = PoolSource::getNextItemType();
        } while (type == ItemType::IsLumi);
      } else {
        seenFirstLumi_ = true;
      }
    }
    return type;
  }
}  // namespace edm
using namespace edm;

DEFINE_FWK_INPUT_SOURCE(OneLumiPoolSource);
