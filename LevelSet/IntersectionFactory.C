/*
 * IntersectionFactory.C
 *
 *  Created on: Feb 12, 2009
 *      Author: Michel Lesoinne
 */

#include <iostream>
#include "IntersectionFactory.h"
#include "Domain.h"

std::map<std::string, IntersectorConstructor *>
IntersectionFactory::allIntersectors = std::map<std::string, IntersectorConstructor *> ();
Communicator *IntersectionFactory::com = 0;

typedef std::map<std::string, IntersectorConstructor *> IntersectorMap;

IntersectorConstructor *lastI = 0;

IntersectorConstructor *
IntersectionFactory::registerClass(std::string name, IntersectorConstructor *o) {
//  std::cout << "Registering Intersector " << name << " at " << o << std::endl;
  o->print();
  lastI = o;
  allIntersectors[name] = o;
  return o;
}

void
IntersectionFactory::parseIntersectionObject(std::string name, ParseTree &data) {
  IntersectorMap::iterator it = allIntersectors.find(name);
  if(it == allIntersectors.end())
    return;
  //return (it == allIntersectors.end()) ? 0 : it->second;
  IntersectorConstructor *ctr = it->second;
  ctr->init(data);
}

DistLevelSetStructure *
IntersectionFactory::getIntersectionObject(std::string name, Domain &domain) {
  IntersectorMap::iterator it = allIntersectors.find(name);
  if(it == allIntersectors.end())
    return 0;
  //return (it == allIntersectors.end()) ? 0 : it->second;
  IntersectorConstructor *ctr = it->second;
  return ctr->getIntersector(*new IntersectProblemData(domain));
}

IntersectProblemData::IntersectProblemData(Domain &d) : domain(d) {}

int IntersectProblemData::getLocalNumSub() {
  return domain.getNumLocSub();
}
