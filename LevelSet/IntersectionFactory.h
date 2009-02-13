#ifndef _INTERSECTION_FACTORY_H_
#define _INTERSECTION_FACTORY_H_

#include <string>
#include <map>

class DistLevelSetStructure;
class ParseTree;
class Communicator;

class IntersectorConstructor {
  public:
    virtual DistLevelSetStructure *getIntersector() = 0;
    virtual void init(ParseTree &dataTree) = 0;
    virtual int print() = 0;
};

/** Factory of Intersection objects
 *
 *  Intersection objects can register with the factory and the
 *  factory informs the parser about names for the available objects.
 *  It also constructs Intersection objects on demand. */
class IntersectionFactory {
    static std::map<std::string, IntersectorConstructor*> allIntersectors;
    static Communicator *com;
  public:

    static IntersectorConstructor *
        registerClass(std::string name, IntersectorConstructor *);
    static DistLevelSetStructure *getIntersectionObject(std::string name, ParseTree &data);
    static Communicator *getCommunicator() { return com; }
    static void setCommunicator(Communicator *c) { com = c; }
};

#endif
