
/*
 * Simple testing of branching lib
 */

#include <iostream>
#include <memory>
#include <vector>

#include "branching.h"

int main(int argc, char const ** argv)
{

  auto ic = std::make_shared<TestCell>();
  std::vector< std::shared_ptr<TestCell> > tst;
  tst.push_back(ic);

  // n cell listener
  auto Nlst = std::make_shared<NCellListener>();

  BProcess<TestCell> bp(1);
  bp.add_listener(Nlst);

  bp.run(100,10);
  
  Nlst->print();

  return 0;
}
