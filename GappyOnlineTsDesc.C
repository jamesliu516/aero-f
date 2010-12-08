#include <DistTimeState.h>
#include <GeoSource.h>
#include <SpaceOperator.h>
#include <Domain.h>
#include <MatVecProd.h>
#include <NewtonSolver.h>

#include <TsOutput.h>

#include <algorithm>
#include <numeric>

//------------------------------------------------------------------------------

template <int dim>
template <typename InputIterator>
RestrictionMapping<dim>::RestrictionMapping(const Domain * domain, InputIterator globalIndexBegin, InputIterator globalIndexEnd) :
  localSubdomainCount_(domain->getNumLocSub()),
  originDistInfo_(domain->getNodeDistInfo()),
  restrictedDistInfo_(originDistInfo().numLocThreads,
                      originDistInfo().numLocSub,
                      originDistInfo().numGlobSub,
                      originDistInfo().locSubToGlobSub,
                      originDistInfo().com),
  sampleNodes_(globalIndexBegin, globalIndexEnd)
{
  originToRestricted_.resize(restrictedDistInfo_.numLocSub);
  restrictedToOrigin_.resize(restrictedDistInfo_.numLocSub);

  for (int iSub = 0; iSub < localSubdomainCount(); ++iSub) {
    const int * nodeTopology = domain->getSubDomain()[iSub]->getNodeMap();

    const int iNodeEnd = originDistInfo().subSize(iSub);
    for (int iNode = 0; iNode < iNodeEnd; ++iNode) {
      const int globalNodeRank = nodeTopology[iNode];
      if (sampleNodes_.find(globalNodeRank) != sampleNodes_.end()) {
        originToRestricted_[iSub].insert(std::make_pair(iNode, restrictedToOrigin_[iSub].size()));
        restrictedToOrigin_[iSub].push_back(iNode);
      }
    }

    restrictedDistInfo_.setLen(iSub, restrictedToOrigin_[iSub].size());
  }

  restrictedDistInfo_.finalize(true);

  for (int iSub = 0; iSub < localSubdomainCount(); ++iSub) {
    const std::vector<int> & subDomainMapping = restrictedToOrigin_[iSub];
    const int nodeCount = subDomainMapping.size();
    for (int iRestNode = 0; iRestNode < nodeCount; ++iRestNode) {
      const int iOrigNode = subDomainMapping[iRestNode]; 
      restrictedDistInfo_.getMasterFlag(iSub)[iRestNode] = originDistInfo().getMasterFlag(iSub)[iOrigNode];
      restrictedDistInfo_.getInvWeight(iSub)[iRestNode] = originDistInfo().getInvWeight(iSub)[iOrigNode];
    }
  }
}

//------------------------------------------------------------------------------

template <int dim>
const DistSVec<double, dim> &
RestrictionMapping<dim>::restriction(const DistSVec<double, dim> & in, DistSVec<double, dim> & out) const {
  assert(&in.info() == &originDistInfo());
  assert(&out.info() == &restrictedDistInfo());

  for (int iSub = 0; iSub < localSubdomainCount(); ++iSub) {
    const SVec<double, dim> & subIn = in(iSub);
    SVec<double, dim> & subOut = out(iSub);

    const std::vector<int> & mapping = restrictedToOrigin_[iSub];
    for (int iRestrict = 0; iRestrict < mapping.size(); ++iRestrict) {
      const int iOrigin = mapping[iRestrict];
      std::copy(subIn[iOrigin], subIn[iOrigin] + dim, subOut[iRestrict]);
    }
  }
  
  return out; 
}

//------------------------------------------------------------------------------

template <int dim>
const DistSVec<double, dim> &
RestrictionMapping<dim>::expansion(const DistSVec<double, dim> & in, DistSVec<double, dim> & out) const {
  assert(&in.info() == &restrictedDistInfo());
  assert(&out.info() == &originDistInfo());
 
  out = 0.0;

  for (int iSub = 0; iSub < localSubdomainCount(); ++iSub) {
    const SVec<double, dim> & subIn = in(iSub);
    SVec<double, dim> & subOut = out(iSub);

    const std::vector<int> & mapping = restrictedToOrigin_[iSub];
    for (int iRestrict = 0; iRestrict < mapping.size(); ++iRestrict) {
      const int iOrigin = mapping[iRestrict];
      std::copy(subIn[iRestrict], subIn[iRestrict] + dim, subOut[iOrigin]);
    }
  }
  
  return out;
}

//------------------------------------------------------------------------------

template <int dim>
double
RestrictionMapping<dim>::dotProduct(const DistSVec<double, dim> & originVec, const DistSVec<double, dim> & restrictedVec) const {
  assert(&originVec.info() == &originDistInfo());
  assert(&restrictedVec.info() == &restrictedDistInfo());

#ifndef MPI_OMP_REDUCTION
  const int globalSubdomainCount = originDistInfo().numGlobSub;
  double * subBuffer = reinterpret_cast<double *>(alloca(sizeof(double) * globalSubdomainCount));

  for (int iSub = 0; iSub < globalSubdomainCount; ++iSub) {
    subBuffer[iSub] = 0.0;
  }
#endif
 
 double result = 0.0;

  if (restrictedDistInfo().masterFlag) {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: result)
#else
#pragma omp parallel for
#endif
    for (int iSub = 0; iSub < localSubdomainCount(); ++iSub) {
      const SVec<double, dim> & subRestrictedVec = restrictedVec(iSub);
      const SVec<double, dim> & subOriginVec = originVec(iSub);

      const bool * subMasterFlag = restrictedDistInfo().getMasterFlag(iSub);

      double subContrib = 0.0;
      const std::vector<int> & mapping = restrictedToOrigin_[iSub];
      for (int iRestrict = 0; iRestrict < mapping.size(); ++iRestrict) {
        if (subMasterFlag[iRestrict]) {
          const int iOrigin = mapping[iRestrict];
          const double nodeContrib = std::inner_product(subRestrictedVec[iRestrict], subRestrictedVec[iRestrict] + dim, subOriginVec[iOrigin], 0.0);
          subContrib += nodeContrib;
        }
      }
#ifdef MPI_OMP_REDUCTION
      result += subContrib;
#else
      const int globalSubdomainRank = restrictedDistInfo().locSubToGlobSub[iSub];
      subBuffer[globalSubdomainRank] += subContrib;
#endif
    }
  
  } else {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: result)
#else
#pragma omp parallel for
#endif
    for (int iSub = 0; iSub < localSubdomainCount(); ++iSub) {
      const SVec<double, dim> & subRestrictedVec = restrictedVec(iSub);
      const SVec<double, dim> & subOriginVec = originVec(iSub);

      double subContrib = 0.0;
      const std::vector<int> & mapping = restrictedToOrigin_[iSub];
      for (int iRestrict = 0; iRestrict < mapping.size(); ++iRestrict) {
        const int iOrigin = mapping[iRestrict];
        const double nodeContrib = std::inner_product(subRestrictedVec[iRestrict], subRestrictedVec[iRestrict] + dim, subOriginVec[iOrigin], 0.0);
        subContrib += nodeContrib;
      }
#ifdef MPI_OMP_REDUCTION
      result += subContrib;
#else
      const int globalSubdomainRank = restrictedDistInfo().locSubToGlobSub[iSub];
      subBuffer[globalSubdomainRank] += subContrib;
#endif
    }

  }

#ifdef MPI_OMP_REDUCTION
  restrictedDistInfo().com->globalSum(1, &result);
#else
  restrictedDistInfo().com->globalSum(globalSubdomainCount, subBuffer);
  result = std::accumulate(subBuffer, subBuffer + globalSubdomainCount, 0.0);
#endif

  return result;
}

//------------------------------------------------------------------------------

template <int dim>
GappyOnlineTsDesc<dim>::GappyOnlineTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  TsDesc<dim>(ioData, geoSource, dom)
{
  /* Finish initialization of TsDesc */
  this->timeState = new DistTimeState<dim>(ioData, this->spaceOp, this->varFcn, this->domain, this->V);
  
  MemoryPool mp;
  this->mmh = this->createMeshMotionHandler(ioData, geoSource, &mp);

  // HACK
  int sampleNodeHack[] = {248, 260, 287};

  restrictionMapping_.reset(new RestrictionMapping<dim>(this->domain, sampleNodeHack + 0, sampleNodeHack + 3));

  // DEBUG
  restrictionMapping_->restrictedDistInfo().print();

  DistSVec<double, dim> originVec(restrictionMapping_->originDistInfo());
  DistSVec<double, dim> restrictedVec(restrictionMapping_->restrictedDistInfo());

  const double fillValue = 2.0;
  const double fillValueSquare = fillValue * fillValue;
  const double expectedDotProduct = (3 * dim) * fillValueSquare;

  originVec = fillValue;
  restrictionMapping_->restriction(originVec, restrictedVec);

  double dotProduct = restrictedVec * restrictedVec;
  assert(dotProduct == expectedDotProduct);
  
  dotProduct = restrictionMapping_->dotProduct(originVec, restrictedVec);
  assert(dotProduct == expectedDotProduct);
 
  restrictionMapping_->expansion(restrictedVec, originVec);
  dotProduct = originVec * originVec;
  assert(dotProduct == expectedDotProduct);

  originVec = 10.0;
  restrictionMapping_->expansion(restrictedVec, originVec);
  dotProduct = originVec * originVec;
  assert(dotProduct == expectedDotProduct);
  
  restrictionMapping_->restriction(originVec, restrictedVec);
  dotProduct = restrictedVec * restrictedVec;
  assert(dotProduct == expectedDotProduct);

  // EXAMPLE OF USE
  typedef VecSet<DistSVec<double, dim> > OnlineMatrix;
  const int reducedBasisSize = 10;
  OnlineMatrix AMatrix(reducedBasisSize, restrictionMapping_->restrictedDistInfo());
  OnlineMatrix BMatrix(reducedBasisSize, restrictionMapping_->restrictedDistInfo());

  DistSVec<double, dim> residual(restrictionMapping_->originDistInfo());
  Vec<double> rhs(reducedBasisSize);
  for (int iReduced = 0; iReduced < reducedBasisSize; ++iReduced) {
    rhs[iReduced] = restrictionMapping_->dotProduct(residual, BMatrix[iReduced]);
  }

  //MatVecProd<dim, dim> * mvp;
  //mvp->evaluate(iterRank, meshPosition, controlVolumes, state, fluxes);
  //mvp->apply(in, out);
}

//------------------------------------------------------------------------------

template <int dim>
int
GappyOnlineTsDesc<dim>::solveNonLinearSystem(DistSVec<double,dim> & U, int iterRank) {
  // TODO
  return 0;
}
