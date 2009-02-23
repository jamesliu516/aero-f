/*
 * LevelSetCriterion.h
 *
 *  Created on: Feb 16, 2009
 *      Author: michel
 */

#ifndef FLUIDTYPECRITERION_H_
#define FLUIDTYPECRITERION_H_
#include <map>

template<class Scalar>
class Vec;

template <class Scalar>
class DistVec;

/** Utility class to give fluid type */
class FluidTypeCriterion {
  public:
    virtual bool isSameFluid(int i, int j) const = 0;
};

class FluidTypeFromLevelSet : public FluidTypeCriterion {
  const Vec<double> &phi;
  public:
    FluidTypeFromLevelSet(const Vec<double> &p) : phi(p) {}
    bool isSameFluid(int i, int j) const;
};

class DistFluidTypeCriterion {
  public:
    virtual FluidTypeCriterion &operator()(int i) const = 0;
};

class DistFluidTypeFromLevelSet : public DistFluidTypeCriterion {
    FluidTypeFromLevelSet **ft;
  public:
    DistFluidTypeFromLevelSet(const DistVec<double> &p);
    FluidTypeCriterion &operator()(int i) const;
};

class DistLevelSetStructure;
class DistFluidTypeFromIntersect : public DistFluidTypeCriterion {
    const DistLevelSetStructure &dLS;
    std::map<int, FluidTypeCriterion &> *ftcMap;
  public:
    DistFluidTypeFromIntersect(const DistLevelSetStructure &dls) : dLS(dls)
      { ftcMap =  new std::map<int, FluidTypeCriterion &>(); }
    FluidTypeCriterion &operator()(int i) const;
};


#endif /* FLUIDTYPECRITERION_H_ */
