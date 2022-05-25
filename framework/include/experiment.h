//
// Created by Thomas Prosser on 25.05.22.
//

#ifndef BEC_LATTICE_KICK_EXPERIMENT_H
#define BEC_LATTICE_KICK_EXPERIMENT_H

namespace bec {

  class Experiment {
  public:
   virtual void execute() = 0;
  };

}

#endif // BEC_LATTICE_KICK_EXPERIMENT_H
