#include <ReinitializeDistanceToWall.h>
#include <LevelSet/LevelSetStructure.h>
#include <Domain.h>
#include <DistVector.h>
#include <FemEquationTerm.h>

//------------------------------------------------------------------------------
//  class ReinitializeDistanceToWall
//
//    Eikonal equation solver used to compute distance to the wall for use in
//    turbulence modeling and other applications for embedded simulations.
//
//  Modified by Sebastian Grimberg, Winter/Spring 2017
//
//------------------------------------------------------------------------------

// Flag determining use of wall distance predictors
#define USE_PREDICTORS

template <int dimLS, int dim>
ReinitializeDistanceToWall<dimLS,dim>::ReinitializeDistanceToWall(IoData &ioData, Domain &domain, SpaceOperator<dim> &spaceOp)
    : iod(ioData), dom(domain), d2wall(domain.getNodeDistInfo()),
    // wall distance
    tag(domain.getNodeDistInfo()),sortedNodes(domain.getNodeDistInfo()),
    isSharedNode(domain.getNodeDistInfo()),
    // distance error estimation
    d2wnm1(domain.getNodeDistInfo()), d2wnm2(domain.getNodeDistInfo()),
    spaceOp(spaceOp),SAsensitiv(domain.getNodeDistInfo()),predictorTag(domain.getNodeDistInfo())
{
  int nSub = dom.getNumLocSub();

  firstCheckedNode = new int[nSub];
  nSortedNodes = new int[nSub];

  // tag shared nodes
  isSharedNode = 0;
  Connectivity *sharedNodes;
  for (int iSub=0; iSub<nSub; iSub++) {
    for (int i = 0; i < dom.getSubDomain()[iSub]->getNumNeighb(); i++) {
      sharedNodes = dom.getSubDomain()[iSub]->getSharedNodes();
      for (int j = 0; j < sharedNodes->num(i); j++)
        isSharedNode(iSub)[(*sharedNodes)[i][j]] = 1;
    }
  }

  // predictors
  predictorTime[0] = -1.0;
  predictorTime[1] = -1.0;
  predictorTime[2] = -1.0;

}

//------------------------------------------------------------------------------

template <int dimLS, int dim>
ReinitializeDistanceToWall<dimLS,dim>::~ReinitializeDistanceToWall()
{

  delete[] firstCheckedNode;
  delete[] nSortedNodes;

}

//------------------------------------------------------------------------------

template <int dimLS, int dim>
int ReinitializeDistanceToWall<dimLS,dim>::ComputeWallFunction(DistLevelSetStructure &LSS,
                                                            DistSVec<double, 3> &X,
                                                            DistGeoState &distGeoState,
                                                            const double t)
{

  // // no updates after t0
  // if (predictorTime[0] >= 0.0) return;

#ifdef USE_PREDICTORS
  int update = UpdatePredictorsCheckTol(LSS, distGeoState, t);
#else
  int update = 3;
#endif

  if (update > 0) {

    // store with previous exact updates
    predictorTime[2] = predictorTime[1];
    predictorTime[1] = predictorTime[0];
    predictorTime[0] = t;
    d2wnm2 = d2wnm1;
    d2wnm1 = d2wall;

    // compute distance in full domain
    if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::NONITERATIVE) {
      PseudoFastMarchingMethod(LSS, X);
    }
    else if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::ITERATIVE) {
      IterativeMethodUpdate(LSS, X);
    }
    else if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::HYBRID) {
      fprintf(stderr,"*** Warning *** Hybrid wall distance is depreciated, using non-iterative instead\n");
      PseudoFastMarchingMethod(LSS, X);
    }
    else {
      fprintf(stderr, " *** Error *** Unknown wall distance method\n");
      exit(1);
    }

#pragma omp parallel for
    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
      for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i) {
        distGeoState(iSub).getDistanceToWall()[i] = d2wall(iSub)[i][0];
      }
    }

    if (update < 3) // reinitialize predictors
      ReinitializePredictors(distGeoState, &LSS, X);
  }

  return update;
}

//------------------------------------------------------------------------------
//  Wall distance reinitialization algorithm methods
//------------------------------------------------------------------------------

// The following is an adaptation of the Fast Marching Method to Embedded
// Turbulent computation. Adam 2012.09
template <int dimLS, int dim>
void ReinitializeDistanceToWall<dimLS,dim>::PseudoFastMarchingMethod(
    DistLevelSetStructure &LSS, DistSVec<double, 3> &X)
{
  d2wall = 1.0e10;
  tag = -1;

  int nSub = dom.getNumLocSub(), isDone = 0, level = 1; // Level 0 (inActive nodes) and 1 (Embedded surface neighbors)
  while (!isDone) { // Tag and update d at every level
    dom.pseudoFastMarchingMethod<1>(tag, X, d2wall, level, 0, sortedNodes,
      nSortedNodes, firstCheckedNode, isSharedNode, &LSS);
    isDone = 1;
    for (int iSub = 0; iSub < nSub; ++iSub) { // I don't think it is a good idea to OMP parallelize this loop. nSub should be small, though!
      if (nSortedNodes[iSub]!=firstCheckedNode[iSub]) {
        isDone = 0;
        break;
      }
    }
    dom.getCommunicator()->globalMin(1, &isDone);
    ++level;
  }

  // dom.getCommunicator()->fprintf(stderr, "There are %d levels\n", level-1);
}

//------------------------------------------------------------------------------

template <int dimLS, int dim>
void ReinitializeDistanceToWall<dimLS,dim>::IterativeMethodUpdate(DistLevelSetStructure &LSS,
                                                              DistSVec<double, 3> &X)
{
  int nSub = dom.getNumLocSub();

  /* NOTE: max level and hybrid method no longer possible with full domain
           iteration (enforce instead using a max distance if desired)
  int max_level = level-1;
  if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::HYBRID &&
      iod.eqs.tc.tm.d2wall.iterativelvl > 1)
    max_level = min(iod.eqs.tc.tm.d2wall.iterativelvl, max_level); */

  double res0, resScale;
  DistSVec<double, 1> resnm1(dom.getNodeDistInfo());
  DistVec<int> tag0(dom.getNodeDistInfo());

  // initialization (levels 0, 1) and first local sweep pass
  d2wall = 1.0e10;
  resnm1 = d2wall;
  tag = -1;
  int level = 1, isDone = 0, it = 1, isConverged = 0, checkConverged = 0;

  dom.pseudoFastMarchingMethodSerial<1>(tag, X, d2wall, level, 0,
    sortedNodes, nSortedNodes, firstCheckedNode, isSharedNode, &LSS);
  dom.pseudoFastMarchingMethodComm<1>(tag, d2wall, sortedNodes, nSortedNodes, it);
  level++;

  while (!isDone) {
    dom.pseudoFastMarchingMethodSerial<1>(tag, X, d2wall, level, 0,
      sortedNodes, nSortedNodes, firstCheckedNode, isSharedNode, &LSS);
    isDone = 1;
    for (int iSub = 0; iSub < nSub; ++iSub) {
      if (nSortedNodes[iSub] != firstCheckedNode[iSub]) {
        isDone = 0;
        break;
      }
    }
    level++;
  }
  it++;

  // for convergence check
  checkConverged = 1;
  for (int iSub = 0; iSub < nSub; ++iSub)
    if (!nSortedNodes[iSub]) { checkConverged = 0; break; } // have not visited processor yet
  if (checkConverged) { // compute scaling residual at first valid iteration for subdomain
    resScale = 0.0;
    for (int iSub = 0; iSub < nSub; ++iSub) {
      for (int i = 0; i < d2wall(iSub).len; i++) {
        resScale += d2wall(iSub)[i][0]*d2wall(iSub)[i][0];
      }
    }
    resScale = 1.0/sqrt(resScale);
  }

  // wipe all tags, set embedded surface neighbors fixed (tag = 0)
  tag0 = -1;
  for (int iSub = 0; iSub < nSub; ++iSub) {
    for (int i = 0; i < d2wall(iSub).len; i++) {
      if (tag(iSub)[i] == 1 || tag(iSub)[i] == 0) tag0(iSub)[i] = 0;
    }
  }

  // iteration loop
  while (!isConverged) {

    // reset tags, nSortedNodes and firstCheckedNode at start of each iteration
    tag = tag0;
    for (int iSub = 0; iSub < nSub; ++iSub) {
      firstCheckedNode[iSub] = 0;
      nSortedNodes[iSub] = 0;
    }

    // share information across processes (min node tag and active list level 1 populate)
    dom.pseudoFastMarchingMethodComm<1>(tag, d2wall, sortedNodes, nSortedNodes, it, &resnm1);

    // sweep subdomains
    resnm1 = d2wall;
    level = 2, isDone = 0;  // active list already populated in comm, tags reset
    while (!isDone) {
      dom.pseudoFastMarchingMethodSerial<1>(tag, X, d2wall, level, 1,
        sortedNodes, nSortedNodes, firstCheckedNode, isSharedNode, &LSS);
      isDone = 1;
      for (int iSub = 0; iSub < nSub; ++iSub) {
        if (nSortedNodes[iSub] != firstCheckedNode[iSub]) {
          isDone = 0;
          break;
        }
      }
      ++level;
    }

    if (checkConverged) { // scaled 2-norm convergence check (local norm)
      res0 = 0.0;
      int nActiveNodes = 0;
      for (int iSub = 0; iSub < nSub; ++iSub) {
        for (int i = 0; i < d2wall(iSub).len; i++) {
          res0 += (d2wall(iSub)[i][0]-resnm1(iSub)[i][0])
                 *(d2wall(iSub)[i][0]-resnm1(iSub)[i][0]);
        }
      }
      res0 = sqrt(res0)*resScale;
      isConverged = (res0 < iod.eqs.tc.tm.d2wall.eps || it > iod.eqs.tc.tm.d2wall.maxIts);
    }
    else {
      checkConverged = 1;
      for (int iSub = 0; iSub < nSub; ++iSub)
        if (!nSortedNodes[iSub]) {checkConverged = 0; break;} // have not visited processor yet
      if (checkConverged) { // compute scaling residual at first valid iteration for subdomain
        resScale = 0.0;
        for (int iSub = 0; iSub < nSub; ++iSub) {
          for (int i = 0; i < d2wall(iSub).len; i++) {
            resScale += d2wall(iSub)[i][0]*d2wall(iSub)[i][0];
          }
        }
        resScale = 1.0/sqrt(resScale);
      }
    }
    dom.getCommunicator()->globalMin(1, &isConverged);
    it++;
  }

  // dom.getCommunicator()->globalSum(1,&res0);
  // res0 /= dom.getCommunicator()->size();
  // dom.getCommunicator()->fprintf(stderr,
  //   "Iterative distance to wall: final residual (CPU avg) = %e, target = %e @ it %d\n\n",
  //   res0, iod.eqs.tc.tm.d2wall.eps, it-1);
}

//------------------------------------------------------------------------------
//  Wall distance error estimator methods
//------------------------------------------------------------------------------

// error code key:  0 : no tolerances exceeded, do not update d2wall
//                  2 : delta tolerance exceeded, full update d
//                  3 : not enough timestep info to initialize predictors, full update d
template <int dimLS, int dim>
int ReinitializeDistanceToWall<dimLS,dim>::UpdatePredictorsCheckTol(
  DistLevelSetStructure &LSS, DistGeoState &distGeoState, const double t)
{
  int update = 0;

  if (predictorTime[1] < 0.0)
    return 3;
  else if (predictorTime[2] < 0.0)
    return 2;

  // compute distance predictions
  double dtnm2 = predictorTime[1]-predictorTime[2];
  double dtnm1 = predictorTime[0]-predictorTime[1];
  double dtn = t-predictorTime[0];

  double dtnodtnm1 = dtn/dtnm1;
  double dtnodtnm2 = dtn/dtnm2;
  double dtnm1odtnm2 = dtnm1/dtnm2;
  double dtnpdtnm1odtnm1pdtnm2 = (dtn+dtnm1)/(dtnm1+dtnm2);

  double d2wp, delta, meandelta = 0.0;  // RMS error

#pragma omp parallel for
  for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
    for (int k = 0; k < d2wall(iSub).len; k++) {
      if (predictorTag(iSub)[k] > 1) {

        // predicted distance
        d2wp = (1.0+dtnodtnm1*(1.0+dtnpdtnm1odtnm1pdtnm2))*d2wall(iSub)[k][0]
          - dtnodtnm1*(dtnodtnm2+dtnm1odtnm2+1.0)*d2wnm1(iSub)[k][0]
          + dtnodtnm2*dtnpdtnm1odtnm1pdtnm2*d2wnm2(iSub)[k][0];

        // local delta definition based on SA source error
        delta = SAsensitiv(iSub)[k]*(d2wp-d2wall(iSub)[k][0]);
        meandelta += delta*delta;
      }
      else if (LSS(iSub).isActive(0.0,k) && predictorTag(iSub)[k] < 1) { // ghost real transition

        distGeoState(iSub).getDistanceToWall()[k] = LSS(iSub).distToInterface(0.0,k);
        predictorTag(iSub)[k] = 1; // modify tag so it is ignored as a predictor node

      }
    }
  }

  // global error metric reduction
  dom.getCommunicator()->globalSum(1,&meandelta);
  meandelta = sqrt(meandelta)/SAsensitivScale;

  if (meandelta > iod.eqs.tc.tm.d2wall.predictoreps)
    update = 2;

  return update;
}

//------------------------------------------------------------------------------

// create new predictors and populate
template <int dimLS, int dim>
void ReinitializeDistanceToWall<dimLS,dim>::ReinitializePredictors(
  DistGeoState &distGeoState, DistLevelSetStructure *LSS, DistSVec<double, 3> &X)
{

  // reset predictors
  SAsensitiv = 0.0;
  SAsensitivScale = 0.0;

  // compute distance sensitivities
  DistSVec<double,dim> *V = spaceOp.getCurrentPrimitiveVector();
  dom.computeSADistSensitivity(spaceOp.getFemEquationTerm(), distGeoState, X,
    *V, SAsensitiv, LSS);

  // retrieve residual for scaling
  DistSVec<double,dim>& RR = *spaceOp.getCurrentSpaceResidual();

#pragma omp parallel for
  for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
    for (int k = 0; k < d2wall(iSub).len; k++) {
      if (d2wall(iSub)[k][0] < 1.0e-15) // currently a ghost node (predictorTag = 0)
        predictorTag(iSub)[k] = 0;
      else if (d2wnm1(iSub)[k][0] < 1.0e-15 || d2wnm2(iSub)[k][0] < 1.0e-15) // node is active but should not be used for predictors
        predictorTag(iSub)[k] = 1;
      else {
        SAsensitivScale += RR(iSub)[k][5]*RR(iSub)[k][5];
        predictorTag(iSub)[k] = 2; // label node as predictor
      }
    }
  }

  // L2 norm of SA residual for predictor nodes (for scaling)
  dom.getCommunicator()->globalSum(1,&SAsensitivScale);
  SAsensitivScale = sqrt(SAsensitivScale);
}

//------------------------------------------------------------------------------
//  Debugging methods
//------------------------------------------------------------------------------

// sjg, 02/2017: instead of computing error for embedded cylinder, compute relative
// error of iterative and noniterative methods
template <int dimLS, int dim>
void ReinitializeDistanceToWall<dimLS,dim>::ComputeExactErrors(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState)
{
  double localError = 0.0;
  int nSub = dom.getNumLocSub();
  int nDofs[nSub];
  double **errors;
  errors = new double *[nSub];

  double localErrorEx = 0.0;
  double **errorsEx;
  errorsEx = new double *[nSub];
  DistSVec<double, 1> d2wall_comp = d2wall;

  // compare to iterative
  int iterativeTemp = iod.eqs.tc.tm.d2wall.maxIts;
  iod.eqs.tc.tm.d2wall.maxIts = 100;
  dom.getCommunicator()->fprintf(stderr,
    "Comparing specified wall distance computation to iterative method with 100 maximum iterations.\n");
  // DistanceToClosestPointOnMovingStructure(LSS, X, distGeoState);
  // GetLevelsFromInterfaceAndMarchForward(LSS, X, distGeoState);
  IterativeMethodUpdate(LSS, X); // new iterative implementation
  iod.eqs.tc.tm.d2wall.maxIts = iterativeTemp;

  // // compare to noniterative
  // dom.getCommunicator()->fprintf(stderr,
  //   "Comparing specified wall distance computation to noniterative method.\n");
  // PseudoFastMarchingMethod(LSS,X);

  DistSVec<double, 1> d2wall_ref = d2wall;

#pragma omp parallel for
  for (int iSub = 0; iSub < nSub; ++iSub)
  {
    errors[iSub] = new double[4];
    errors[iSub][0] = 0.0;          // mean
    errors[iSub][1] = 0.0;          // RMS
    errors[iSub][2] = 0.0;          // max
    errors[iSub][3] = 0.0;          // distance of max
    errorsEx[iSub] = new double[4];
    errorsEx[iSub][0] = 0.0;
    errorsEx[iSub][1] = 0.0;
    errorsEx[iSub][2] = 0.0;
    errorsEx[iSub][3] = 0.0;

    nDofs[iSub] = 0;
    for (int i = 0; i < d2wall(iSub).size(); ++i)
    {
      if (LSS(iSub).isActive(0.0, i))
      {
        nDofs[iSub]++;
        localError = fabs(d2wall_ref(iSub)[i][0] - d2wall_comp(iSub)[i][0]);
        errors[iSub][0] += localError;
        errors[iSub][1] += localError * localError;
        if (localError > errors[iSub][2])
        {
          errors[iSub][2] = localError;
          errors[iSub][3] = d2wall_ref(iSub)[i][0];
        }
        localErrorEx = fabs(d2wall_ref(iSub)[i][0] - d2wall_comp(iSub)[i][0]) / d2wall_ref(iSub)[i][0];
        errorsEx[iSub][0] += localErrorEx;
        errorsEx[iSub][1] += localErrorEx * localErrorEx;
        if (localErrorEx > errorsEx[iSub][2])
        {
          errorsEx[iSub][2] = localErrorEx;
          errorsEx[iSub][3] = d2wall_ref(iSub)[i][0];
        }
      }
    }
  }
  for (int iSub = 1; iSub < nSub; ++iSub)
  {
    nDofs[0] += nDofs[iSub];
    errors[0][0] += errors[iSub][0];
    errors[0][1] += errors[iSub][1];
    if (errors[iSub][2] > errors[0][2])
    {
      errors[0][2] = errors[iSub][2];
      errors[0][3] = errors[iSub][3];
    }
    errorsEx[0][0] += errorsEx[iSub][0];
    errorsEx[0][1] += errorsEx[iSub][1];
    if (errorsEx[iSub][2] > errorsEx[0][2])
    {
      errorsEx[0][2] = errorsEx[iSub][2];
      errorsEx[0][3] = errorsEx[iSub][3];
    }
  }

  // Communicate across all processes to find global sums/max
  dom.getCommunicator()->globalSum(1, nDofs);
  dom.getCommunicator()->globalSum(2, errors[0]);
  dom.getCommunicator()->globalSum(2, errorsEx[0]);

  // NOTE: incorrect, MPI_2DOUBLE_PRECISION is only defined for FORTRAN codes, use MPI_DOUBLE_INT for C++
  MPI_Comm comm = dom.getCommunicator()->getMPIComm();
  MPI_Allreduce(&errors[0][2], &errors[0][2], 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, comm);
  MPI_Allreduce(&errorsEx[0][2], &errorsEx[0][2], 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, comm);

  errors[0][0] /= nDofs[0];
  errors[0][1] /= nDofs[0];
  errors[0][1] = sqrt(errors[0][1]);
  errorsEx[0][0] /= nDofs[0];
  errorsEx[0][1] /= nDofs[0];
  errorsEx[0][1] = sqrt(errorsEx[0][1]);

  dom.getCommunicator()->fprintf(stderr, "Absolute d2wall Error: %12.8e, %12.8e, %12.8e at %12.8e\n",
    errors[0][0], errors[0][1], errors[0][2], errors[0][3]);
  dom.getCommunicator()->fprintf(stderr, "Relative d2wall Error: %12.8e, %12.8e, %12.8e at %12.8e\n\n",
    errorsEx[0][0], errorsEx[0][1], errorsEx[0][2], errorsEx[0][3]);

  for (int iSub = 0; iSub < nSub; iSub++)
    delete[] errors[iSub];
  delete[] errors;
  for (int iSub = 0; iSub < nSub; iSub++)
    delete[] errorsEx[iSub];
  delete[] errorsEx;
  return;
}

//------------------------------------------------------------------------------

template <int dimLS, int dim>
void ReinitializeDistanceToWall<dimLS,dim>::PrescribedValues(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState)
{
  double mind = 1e10, maxd = -1e10;
#pragma omp parallel for
  for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub)
  {
    for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i)
    {
      d2wall(iSub)[i][0] = distGeoState(iSub).getDistanceToWall()[i];
      mind = min(mind, d2wall(iSub)[i][0]);
      maxd = max(maxd, d2wall(iSub)[i][0]);
    }
  }
  dom.getCommunicator()->globalMin(1, &mind);
  dom.getCommunicator()->globalMax(1, &maxd);
  dom.getCommunicator()->fprintf(stderr, "Min: %e\t\tMax: %e\n", mind, maxd);
  return;
}

//------------------------------------------------------------------------------

template class ReinitializeDistanceToWall<1,5>;
template class ReinitializeDistanceToWall<1,6>;
template class ReinitializeDistanceToWall<1,7>;
template class ReinitializeDistanceToWall<2,5>;
template class ReinitializeDistanceToWall<2,6>;
template class ReinitializeDistanceToWall<2,7>;