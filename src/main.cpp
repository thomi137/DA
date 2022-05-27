//
// Created by Thomas Prosser on 25.05.22.
//

#include <algorithm>
#include <iostream>
#include <iterator>

#include "options.hpp"
#include "bec_groundstate.h"

int main(int argc, char* argv[]) {
  bec::BecOptions options = bec::BecOptions(argc, argv);
  // bec::BecGroundstate exp;
  bec::BecGroundstateImTime exp;
  exp.execute();

}