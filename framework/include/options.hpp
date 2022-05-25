//
// Created by Thomas Prosser on 25.05.22.
//

#ifndef BEC_LATTICE_KICK_OPTIONS_HPP
#define BEC_LATTICE_KICK_OPTIONS_HPP

#include <iostream>
#include <boost/program_options.hpp>

#include "bec.h"
namespace po = boost::program_options;

namespace bec {

  class BecOptions {
  public:
    BecOptions(int argc, char** argv): argc_(argc), argv_(argv) {
      init();
    }

    void init() {
      try {
        po::options_description desc("Options for controlling a kicked BEC");
        desc.add_options()
            ("help,h", "Produce a help screen")
            ("kick_strength,K", po::value<bec_t>()->default_value(0.0), "Strength of kick")
            ("g-coupling,g", po::value<bec_t>()->default_value(0.0), "Interaction coupling for BEC")
            ("lattice,l", po::value<bool>()->default_value(false), "lattice")
            ("potential,p", po::value<bool>()->default_value(false), "potential")
            ;
        po::store(po::parse_command_line(argc_, argv_, desc), vm_);
        po::notify(vm_);

        if (vm_.count("help")) {
          std::cout << "help set" << std::endl;
        }

        if(vm_.count("kick_strength")) {
          printf("Kick strength: %1f\n", vm_["kick_strength"].as<bec_t>());
        }

      } catch (std::exception& e){
        std::cerr << e.what() << std::endl;
      }
    }

   po::variables_map get_variables_map() {return vm_;}

  private:
    int argc_;
    char** argv_;
    po::variables_map vm_;
  };

}

#endif // BEC_LATTICE_KICK_OPTIONS_HPP
