#ifndef _GAPPY_ONLINE_TS_DESC_H_
#define _GAPPY_ONLINE_TS_DESC_H_

#include <IoData.h>
#include <TsDesc.h>
#include <KspPrec.h>

struct DistInfo;

class GeoSource;
class Domain;
class Communicator;

template <class Scalar, int dim> class DistSVec;
template <int dim, int neq> class MatVecProdFD;

#include <vector>
#include <set>
#include <map>
#include <memory>

//------------------------------------------------------------------------------

template <int dim>
class RestrictionMapping {
public:
  int localSubdomainCount() const { return localSubdomainCount_; }
  const DistInfo & originDistInfo() const { return originDistInfo_; }
  const DistInfo & restrictedDistInfo() const { return restrictedDistInfo_; }

  const DistSVec<double, dim> & restriction(const DistSVec<double, dim> &, DistSVec<double, dim> &) const;
  const DistSVec<double, dim> & expansion(const DistSVec<double, dim> &, DistSVec<double, dim> &) const;
 
  double dotProduct(const DistSVec<double, dim> & originVec, const DistSVec<double, dim> & restrictedVec) const;

  template <typename InputIterator>
  RestrictionMapping(const Domain * domain, InputIterator globalIndexBegin, InputIterator globalIndexEnd);

private:
  typedef std::map<int, int> NumberingMap;

  int localSubdomainCount_;
  const DistInfo & originDistInfo_;
  DistInfo restrictedDistInfo_;
  
  std::set<int> sampleNodes_;
  std::vector<NumberingMap> originToRestricted_;
  std::vector<std::vector<int> > restrictedToOrigin_;

  // Disallow copy and assignment
  RestrictionMapping(const RestrictionMapping &); // = delete;
  const RestrictionMapping & operator=(const RestrictionMapping &); // = delete;
};

//------------------------------------------------------------------------------

template <int dim>
class GappyOnlineTsDesc : public TsDesc<dim> {
public:
  virtual int solveNonLinearSystem(DistSVec<double,dim> & U, int iterRank); // overriden

  GappyOnlineTsDesc(IoData &, GeoSource &, Domain *);

private:
  std::auto_ptr<RestrictionMapping<dim> > restrictionMapping_;
  
  // Disallow copy and assignment
  GappyOnlineTsDesc(const GappyOnlineTsDesc &); // = delete;
  const GappyOnlineTsDesc & operator=(const GappyOnlineTsDesc &); // = delete;
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <GappyOnlineTsDesc.C>
#endif

#endif
