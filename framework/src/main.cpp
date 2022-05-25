//
// Created by Thomas Prosser on 25.05.22.
//

#include <algorithm>
#include <iostream>
#include <iterator>

#include "../include/options.hpp"
#include "../include/bec_groundstate/bec_groundstate.h"

int main(int argc, char* argv[]) {
  bec::BecOptions options = bec::BecOptions(argc, argv);
  bec::BecGroundstate exp;
  exp.execute();

}