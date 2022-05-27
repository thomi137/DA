//
// Created by Thomas Prosser on 25.05.22.
//

#ifndef BEC_LATTICE_KICK_BEC_GROUNDSTATE_H
#define BEC_LATTICE_KICK_BEC_GROUNDSTATE_H

#include "experiment.h"

namespace bec {

  class BecGroundstate : public Experiment {
  public:
    ~BecGroundstate() override = default;
    void execute() override;
  };

  class BecGroundstateImTime : public Experiment {
  public:
    ~BecGroundstateImTime() override = default;
    void execute() override;
  };


} // namespace bec

#endif // BEC_LATTICE_KICK_BEC_GROUNDSTATE_H
