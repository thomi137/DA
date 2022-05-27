//
// Created by Thomas Prosser on 25.05.22.
//

#ifndef BEC_LATTICE_KICK_EXPERIMENT_H
#define BEC_LATTICE_KICK_EXPERIMENT_H

namespace bec {

  class Experiment {
  public:
    virtual ~Experiment() = default;
    virtual void execute() = 0;
  };

} // namespace bec

#endif // BEC_LATTICE_KICK_EXPERIMENT_H
