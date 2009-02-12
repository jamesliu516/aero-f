/*
 * IntersectionFactory.C
 *
 *  Created on: Feb 12, 2009
 *      Author: Michel Lesoinne
 */

#include <iostream>
#include "IntersectionFactory.h"

std::map<std::string, IntersectorConstructor *>
IntersectionFactory::allIntersectors = std::map<std::string, IntersectorConstructor *> ();

typedef std::map<std::string, IntersectorConstructor *> IntersectorMap;

IntersectorConstructor *lastI = 0;

IntersectorConstructor *
IntersectionFactory::registerClass(std::string name, IntersectorConstructor *o) {
  std::cout << "Registering Intersector " << name << " at " << o << std::endl;
  o->print();
  lastI = o;
  allIntersectors[name] = o;
  return o;
}

DistLevelSetStructure *
IntersectionFactory::getIntersectionObject(std::string name, ParseTree &data) {
  IntersectorMap::iterator it = allIntersectors.find(name);
  if(it == allIntersectors.end())
    return 0;
  //return (it == allIntersectors.end()) ? 0 : it->second;
  IntersectorConstructor *ctr = it->second;
  std::cout << "Got the intersector at " << ctr << std::endl;
  lastI->print();
  ctr->print();
  ctr->init(data);
  return ctr->getIntersector();
}
