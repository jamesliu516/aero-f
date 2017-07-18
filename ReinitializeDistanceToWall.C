#include <ReinitializeDistanceToWall.h>
#include <LevelSet/LevelSetStructure.h>
#include <Domain.h>
#include <DistVector.h>
#include <FemEquationTerm.h>

//------------------------------------------------------------------------------
//  ReinitializeDistanceToWall class
//
//    Eikonal equation solver used to compute distance to the wall for use in
//    turbulence modeling and other applications for embedded simulations.
//
//  Modified by Sebastian Grimberg, Winter/Spring 2017
//------------------------------------------------------------------------------

// // flags determining use of wall distance predictors
//   #define USE_PREDICTORS
//   #define PREDICTOR_DEBUG

template <int dimLS, int dim>
ReinitializeDistanceToWall<dimLS,dim>::ReinitializeDistanceToWall(IoData &ioData, Domain &domain, SpaceOperator<dim> &spaceOp)
    : iod(ioData), dom(domain), d2wall(domain.getNodeDistInfo()),
    // default wall distance
    tag(domain.getNodeDistInfo()),sortedNodes(domain.getNodeDistInfo()),
    // FEM wall distance
    nodeTag(domain.getNodeDistInfo()),isSharedNode(domain.getNodeDistInfo()),
    unsortedTag(domain.getNodeDistInfo()),unsortedNodes(domain.getNodeDistInfo()),
    // predictors
    d2wnm1(domain.getNodeDistInfo()), d2wnm2(domain.getNodeDistInfo()), countReinits(0), nPredTot(0),
    spaceOp(spaceOp),SAsensi(domain.getNodeDistInfo())
{
  int nSub = dom.getNumLocSub();

  // default wall distance
  // nActiveNodes = new int[nSub];
  firstCheckedNode = new int[nSub];
  nSortedNodes = new int[nSub];

  // FEM wall distance
  nSortedElems = new int[nSub];
  firstCheckedElem = new int[nSub];
  nUnsortedNodes = new int[nSub];

  activeElemList = new int *[nSub];
  elemTag = new int *[nSub];
  knownNodes = new int *[nSub];
  for (int iSub=0; iSub<nSub; iSub++) {
    activeElemList[iSub] = new int[dom.getSubDomain()[iSub]->numElems()];
    elemTag[iSub] = new int[dom.getSubDomain()[iSub]->numElems()];
    knownNodes[iSub] = new int[dom.getSubDomain()[iSub]->numElems()];
  }

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
  SAsensi = 0.0;

  tolmax = 5.0e-2, tolmin = 1.0e-8; // predictor tolerances, to be read from parser eventually
}

//------------------------------------------------------------------------------

template <int dimLS, int dim>
ReinitializeDistanceToWall<dimLS,dim>::~ReinitializeDistanceToWall()
{

#ifdef PREDICTOR_DEBUG
  dom.getCommunicator()->fprintf(stderr,"\nWall distance computer: total full domain reinitializations = %d\n\n",
    countReinits);
#endif

  // default wall distance
  // delete[] nActiveNodes;
  delete[] firstCheckedNode;
  delete[] nSortedNodes;

  // FEM wall distance
  delete[] nSortedElems;
  delete[] firstCheckedElem;
  delete[] nUnsortedNodes;

  for (int iSub=0; iSub<dom.getNumLocSub(); iSub++) {
    delete[] activeElemList[iSub];
    delete[] elemTag[iSub];
    delete[] knownNodes[iSub];
  }
  delete[] activeElemList;
  delete[] elemTag;
  delete[] knownNodes;
}

//------------------------------------------------------------------------------

template <int dimLS, int dim>
void ReinitializeDistanceToWall<dimLS,dim>::ComputeWallFunction(DistLevelSetStructure &LSS,
                                                            DistSVec<double, 3> &X,
                                                            DistGeoState &distGeoState,
                                                            const double t)
{

#ifdef USE_PREDICTORS
  // update predictors and check tolerances, receiving global maximum
  int update = UpdatePredictorsCheckTol(LSS, distGeoState, t);
#else
  int update = 3;
#endif

  if (update >= 1) {
    // compute distance in full domain

#ifdef PREDICTOR_DEBUG
    countReinits++;
    dom.getCommunicator()->fprintf(stderr, "Performing a full domain wall distance update!\n",update);
#endif

#ifdef USE_PREDICTORS
    // store with previous exact updates
    predictorTime[2] = predictorTime[1];
    predictorTime[1] = predictorTime[0];
    predictorTime[0] = t;
    d2wnm2 = d2wnm1;
    d2wnm1 = d2wall;
#endif

    // update d2wall exactly
    if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::NONITERATIVE) {
      PseudoFastMarchingMethod(LSS, X, 0); // default
      // PseudoFastMarchingMethodFEM(LSS, X); // new FEM
    }
    else if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::ITERATIVE) {
      // InitializeWallFunction(LSS, X, distGeoState, t);
      // GetLevelsFromInterfaceAndMarchForward(LSS, X, distGeoState);

      // new iterative implementation
      IterativeMethodUpdate(LSS, X);
    }
    else if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::HYBRID) {

      fprintf(stderr,"*** Warning *** Hybrid wall distance is depreciated, calling non-iterative instead\n");
      PseudoFastMarchingMethod(LSS, X, 0);

      // int iterativeLevel = 0;

      // if (iod.eqs.tc.tm.d2wall.iterativelvl > 1) {
      //   // InitializeWallFunction(LSS, X, distGeoState, t);
      //   // GetLevelsFromInterfaceAndMarchForward(LSS, X, distGeoState);

      //   IterativeMethodUpdate(LSS, X);
      //   iterativeLevel = iod.eqs.tc.tm.d2wall.iterativelvl;
      // }
      // PseudoFastMarchingMethod(LSS, X, iterativeLevel);
    }
    else {
      fprintf(stderr, " *** Error *** Unknown wall distance method\n");
      exit(1);
    }

#ifdef USE_PREDICTORS
    // reinitialize predictors
    if (update < 3)
      ReinitializePredictors(update, t, &LSS, X);
#endif

    // ComputePercentChange(LSS, distGeoState);

#pragma omp parallel for
    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
      for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i) {
        distGeoState(iSub).getDistanceToWall()[i] = d2wall(iSub)[i][0];
      }
    }
  }

#ifdef PREDICTOR_DEBUG
  else {
    dom.getCommunicator()->fprintf(stderr, "SUCCESS: Skipped wall distance calculation!\n\n",update);
  }
  // dom.getCommunicator()->fprintf(stderr, "Times: tn = %e, tnm1 = %e, tnm2 = %e\n\n",predictorTime[0],predictorTime[1],predictorTime[2]);
#endif

  // PrintIntersectedValues(&LSS, X, distGeoState);
  // ComputeExactErrors(LSS, X, distGeoState);

  return;
}

//------------------------------------------------------------------------------

// sjg, 06/2017: FEM based wall distance via a marching method
template <int dimLS, int dim>
void ReinitializeDistanceToWall<dimLS,dim>::PseudoFastMarchingMethodFEM(
    DistLevelSetStructure &LSS, DistSVec<double, 3> &X)
{
  d2wall = 1.0e10;
  nodeTag = -1;
  unsortedTag = -1;

  int nSub = dom.getNumLocSub();
  for (int iSub = 0; iSub < nSub; iSub++) {
    std::memset(elemTag[iSub], -1, sizeof(int)*dom.getSubDomain()[iSub]->numElems());
    std::memset(knownNodes[iSub], 0, sizeof(int)*dom.getSubDomain()[iSub]->numElems());
  }

  int isDone = 0, level = 1; // Level 0 (inactive nodes) and 1 (embedded surface neighbors)
  while (!isDone) { // Tag and update d at every level
    dom.pseudoFastMarchingMethodFEM<1>(X, d2wall, nodeTag, level, elemTag, activeElemList,
      knownNodes, nSortedNodes, nSortedElems, firstCheckedElem,
      unsortedTag, unsortedNodes, nUnsortedNodes,
      isSharedNode, &LSS);

    isDone = 1;
    for (int iSub = 0; iSub < nSub; ++iSub) {
      if (nSortedElems[iSub] != firstCheckedElem[iSub]) { // while active list is not empty
        isDone = 0;
        break;
      }
    }
    dom.getCommunicator()->globalMin(1, &isDone);
    ++level;

    // // debug
    // if (level >= 200) break;
  }

  // // debug
  // fprintf(stderr,"Before finalization (level %d): "
  //   "Total sorted nodes = %d out of %d, unsorted nodes = %d (nSortedElems = %d "
  //   "out of %d)\n",level-1,nSortedNodes[0],d2wall(0).len,nUnsortedNodes[0],
  //   nSortedElems[0],dom.getSubDomain()[0]->numElems());

  // int sortedel = 0, sortedel2 = 0;
  // for (int iSub=0; iSub < nSub; iSub++) {
  //   for (int i = 0; i < dom.getSubDomain()[iSub]->numElems(); i++) {
  //     if (elemTag[iSub][i] >= 0) sortedel++;
  //     if (knownNodes[iSub][i] == 4) sortedel2++;
  //   }
  // }
  // fprintf(stderr,"Sum tag>-1 = %d, sum knownNodes=4 = %d\n",sortedel,sortedel2);

  // // finalize missed nodes
  // dom.pseudoFastMarchingMethodFinalize<1>(X, d2wall, knownNodes, nSortedNodes,
  //   unsortedNodes, nUnsortedNodes, nodeTag, unsortedTag,
  //   isSharedNode, &LSS);

  // // debug (check completeness)
  // isDone = 1;
  // for (int iSub = 0; iSub < nSub; ++iSub) {
  //   if (nSortedNodes[iSub] != d2wall(iSub).len) {
  //     isDone = 0;
  //     break;
  //   }
  // }
  // if (isDone != 1) {
  //   fprintf(stderr,"PROBLEM: At end of update (level %d): "
  //     "Total sorted nodes = %d out of %d (nSortedElems = %d "
  //     "out of %d)\n",level-1,nSortedNodes[0],d2wall(0).len,nSortedElems[0],
  //     dom.getSubDomain()[0]->numElems());
  //   exit(-1);
  // }

  // dom.getCommunicator()->barrier();
  // exit(-1);

  // dom.getCommunicator()->fprintf(stderr, "There are %d levels\n", level-1);
}

//------------------------------------------------------------------------------

// The following is an adaptation of the Fast Marching Method to Embedded Turbulent computation.
// Adam 2012.09
template <int dimLS, int dim>
void ReinitializeDistanceToWall<dimLS,dim>::PseudoFastMarchingMethod(
    DistLevelSetStructure &LSS, DistSVec<double, 3> &X, int iterativeLevel)
{
  // sortedNodes = -1;
  // if (iterativeLevel == 0)
  // {
  d2wall = 1.0e10;
  tag = -1;
  // }

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

  // double res0, res, resnm1, resGlob, resnm1Glob;
  // testing
  double res, res0, resScaleTmp, resScale = -1.0; // res is not needed
  DistSVec<double, 1> resnm1(dom.getNodeDistInfo());

  // initialization (levels 0, 1) and first local sweep pass
  d2wall = 1.0e10;
  tag = -1;
  res = 0.0;

  int level = 1, it = 1, isDone = 0, isConverged = 0;
  res += dom.pseudoFastMarchingMethodSerial<1>(tag, X, d2wall, level, 0,
    sortedNodes, nSortedNodes, firstCheckedNode, isSharedNode, &LSS);
  dom.pseudoFastMarchingMethodComm<1>(tag, d2wall, sortedNodes, nSortedNodes, it);
  level++;

  while (!isDone) {
    res += dom.pseudoFastMarchingMethodSerial<1>(tag, X, d2wall, level, 0,
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

  // resGlob = res;
  // dom.getCommunicator()->globalSum(1, &resGlob);
  // resGlob = sqrt(resGlob);

  // wipe all tags, set embedded surface neighbors fixed (tag = 0)
  // ONLY do if nCPU > 1 (otherwise will mess up multiple passes for 1 subdomain)
  if (dom.getNumLocSub() > 1 || dom.getCommunicator()->size() > 1) {
    for (int iSub = 0; iSub < nSub; ++iSub) {
      for (int i = 0; i < d2wall(iSub).len; i++) {
        if (tag(iSub)[i] == 1) tag(iSub)[i] = 0;
        else if (tag(iSub)[i] != 0) tag(iSub)[i] = -1;
      }
    }
  }

  /* NOTE: max level and hybrid method no longer possible with full domain
           iteration (enforce instead using a max distance if desired)
  int max_level = level-1;
  if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::HYBRID &&
      iod.eqs.tc.tm.d2wall.iterativelvl > 1)
    max_level = min(iod.eqs.tc.tm.d2wall.iterativelvl, max_level); */

  // iteration loop
  const double eps = 1.0e-10;
  while (!isConverged) {

    // resnm1Glob = resGlob;
    // resnm1 = res;
    // res = 0.0;
    // testing
    resnm1 = d2wall;

    // share information across processes (min node tag)
    dom.pseudoFastMarchingMethodComm<1>(tag, d2wall, sortedNodes, nSortedNodes, it);

    // sweep subdomains (reset all tags > 1 to -1 on first level = 1 call)
    level = 1, isDone = 0;
    while (!isDone) {
      res += dom.pseudoFastMarchingMethodSerial<1>(tag, X, d2wall, level, 1,
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

    // // check convergence
    // if (res < eps) res = resnm1;
    // if (res < eps) isConverged = 0;
    // else isConverged = 1;
    // dom.getCommunicator()->globalMin(1, &isConverged);

    // if (isConverged) {
    //   resGlob = res;
    //   dom.getCommunicator()->globalSum(1, &resGlob);
    //   resGlob = sqrt(resGlob);
    //   res0 = fabs(resGlob-resnm1Glob)/(resGlob+resnm1Glob);
    //   isConverged = (res0 < iod.eqs.tc.tm.d2wall.eps || it > iod.eqs.tc.tm.d2wall.maxIts);
    // }

    // alteratively, 2 norm convergence check with additional loop over nodes
    resScaleTmp = 0.0;
    res0 = 0.0;
    isConverged = 1;
    for (int iSub = 0; iSub < nSub; ++iSub) {
      for (int i = 0; i < d2wall(iSub).len; i++) {
        if (d2wall(iSub)[i][0] == 1.0e10) {  // not converged: have not visited all subdomains yet
          isConverged = 0;
          goto breakpt;
        }
        resScaleTmp += d2wall(iSub)[i][0]*d2wall(iSub)[i][0];
        res0     += (d2wall(iSub)[i][0]-resnm1(iSub)[i][0])
                      *(d2wall(iSub)[i][0]-resnm1(iSub)[i][0]);
      }
    }
    breakpt:
    dom.getCommunicator()->globalMin(1, &isConverged);

    if (isConverged) {
      if (resScale <= eps) {
        // for (int iSub = 0; iSub < nSub; ++iSub) {
        //   for (int i = 0; i < d2wall(iSub).len; i++) {
        //     resScale += d2wall(iSub)[i][0]*d2wall(iSub)[i][0];
        //   }
        // }
        resScale = resScaleTmp;
        dom.getCommunicator()->globalSum(1, &resScale);
        resScale = sqrt(resScale);
      }

      dom.getCommunicator()->globalSum(1, &res0);
      res0 = sqrt(res0)/resScale;
      isConverged = (res0 < iod.eqs.tc.tm.d2wall.eps || it > iod.eqs.tc.tm.d2wall.maxIts);
    }

    dom.getCommunicator()->fprintf(stderr,"Iteration %d of wall distance (res = %e, isConverged = %d)\n",it,res0,isConverged);

    it++;
  }



  // // int *checkConverged = new int[nSub+1];
  // // for (int iSub = 0; iSub <= nSub; ++iSub) checkConverged[iSub] = 0;

  // while (!isConverged) {
  //   isConverged = 1;
  //   resnm1 = res;

  //   dom.getCommunicator()->fprintf(stderr,"\nIteration %d of wall distance (res = %e)\n\n",it,res0);

  //   // // total wipe of tags after first pass from wall ----------- TESTING!
  //   // if (it == 2) {
  //   //   for (int iSub = 0; iSub < nSub; ++iSub) {
  //   //     for (int i = 0; i < d2wall(iSub).len; i++) {
  //   //       if (tag(iSub)[i] != 0) tag(iSub)[i] = -1;
  //   //     }
  //   //     // nSortedNodes[iSub] = 0;
  //   //     // firstCheckedNode[iSub] = 0;
  //   //   }
  //   // }

  //   // share information across processes
  //   dom.pseudoFastMarchingMethodComm<1>(tag, d2wall, sortedNodes, nSortedNodes, level);

  //   fprintf(stderr,"After comm, nSortedNodes = %d (firstCheckedNode = %d on cpu %d)\n",
  //     nSortedNodes[0],firstCheckedNode[0], dom.getCommunicator()->cpuNum());

  //   // propagate information outwards from minimum level
  //   level = 1, isDone = 0;
  //   // level = 2, isDone = 0;
  //   res = 0.0;

  //   while (!isDone) {
  //     res += dom.pseudoFastMarchingMethodSerial<1>(tag, X, d2wall, level, 1,
  //       sortedNodes, nSortedNodes, firstCheckedNode, isSharedNode, &LSS);
  //     isDone = 1;
  //     for (int iSub = 0; iSub < nSub; ++iSub) {
  //       if (nSortedNodes[iSub]==0)
  //         isConverged = 0;
  //       else if (nSortedNodes[iSub]!=firstCheckedNode[iSub]) {
  //       // if (nSortedNodes[iSub]!=firstCheckedNode[iSub]) {
  //         isDone = 0;
  //         break;
  //       }
  //     }
  //     ++level;
  //   }

  //   // debug
  //   if (it==2) return;

  //   // // convergence  check
  //   // if (checkConverged[nSub] != nSub) {
  //   //   for (int iSub = 0; iSub < nSub; ++iSub) {
  //   //     if (nSortedNodes[iSub] > 0 && !checkConverged[iSub]) {
  //   //       checkConverged[iSub] = 1;
  //   //       checkConverged[nSub]++;
  //   //     }
  //   //   }
  //   // }
  //   // else isConverged = 1;
  //   dom.getCommunicator()->globalMin(1, &isConverged);

  //   dom.getCommunicator()->globalSum(1, &res);
  //   res = sqrt(res);

  //   if (isConverged) { // update residual
  //     res0 = fabs(res-resnm1)/(res+resnm1);
  //     isConverged = (res0 < iod.eqs.tc.tm.d2wall.eps || it > iod.eqs.tc.tm.d2wall.maxIts);

  //     // dom.getCommunicator()->fprintf(stderr, "Residual = %e @ iteration %d (isConverged = %d)\n", res[0], it, isConverged);
  //   }
  //   it++;
  // }

  // delete[] checkConverged;

  dom.getCommunicator()->fprintf(stderr,
    "Iterative distance to wall computation: final residual = %e, target = %e @ it %d\n\n",
    res0, iod.eqs.tc.tm.d2wall.eps, it-1);
}

//------------------------------------------------------------------------------

// template <int dimLS, int dim>
// void ReinitializeDistanceToWall<dimLS,dim>::InitializeWallFunction(DistLevelSetStructure &LSS,
//                                                               DistSVec<double, 3> &X,
//                                                               DistGeoState &distGeoState,
//                                                               int *nPredLoc, double t)
// {
//   sortedNodes = -1;
//   tag = -1;

//   // Fill with initial guess
// #if 1
//   d2wall = 1e10;
// #else
//   if (t>0.0) {
//   #pragma omp parallel for
//     for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub)
//     {
//       for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i)
//         d2wall(iSub)[i][0] = distGeoState(iSub).getDistanceToWall()[i]>0.0?distGeoState(iSub).getDistanceToWall()[i]:1e10;
//     }
//   }
//   else d2wall = 1e10;
// #endif

//   int level = 0; // Level 0 (inActive nodes) and 1 (Embedded surface neighbors)
//   while (level < 2)
//   { // Tag levels 0 (inactive) and 1 (embedded surface neighbors) for FIM initialization
//     dom.pseudoFastMarchingMethod<1>(tag, X, d2wall, level, 0, sortedNodes, nSortedNodes, firstCheckedNode, nPredLoc, &LSS);
//     ++level;
//   }

//   return;
// }

//------------------------------------------------------------------------------

// template <int dimLS, int dim>
// void ReinitializeDistanceToWall<dimLS,dim>::GetLevelsFromInterfaceAndMarchForward(DistLevelSetStructure &LSS,
//                                                                               DistSVec<double, 3> &X,
//                                                                               DistGeoState &distGeoState)
// {
//   int max_level = 1;

//   // int min_level = 0;
//   // int level = 1;
//   int min_level = -1;
//   int level = 2;  // levels <= 1 tagged in initialization already

//   dummyPhi = 1.0;

//   int counter = 0;
//   // Tag every level
//   // while (min_level <= 0)
//   while (min_level < 0)
//   {
//     dom.TagInterfaceNodes(0, tag, dummyPhi, level, &LSS);
//     min_level = 1;
//     counter++;

//     for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub)
//     {
//       // for (int i = 0; i < done(iSub).len; ++i)
//       for (int i = 0; i < d2wall(iSub).len; ++i)
//       {
//         min_level = min(min_level, tag[iSub][i]);
//         if (min_level < 0)
//           goto exitlbl;
//         max_level = max(max_level, tag[iSub][i]);
//       }
//     }
//     exitlbl:
//     dom.getCommunicator()->globalMin(1, &min_level);
//     // dom.getCommunicator()->fprintf(stderr, "min_level = %d (count = %d)\n",min_level,counter);
//     ++level;
//   }
//   dom.getCommunicator()->globalMax(1, &max_level);
//   // dom.getCommunicator()->fprintf(stderr, "min_level = %d, max_level = %d\n",min_level,max_level);

//   if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::HYBRID &&
//       iod.eqs.tc.tm.d2wall.iterativelvl > 1)
//     max_level = min(iod.eqs.tc.tm.d2wall.iterativelvl, max_level);

//   // Propagate information outwards
//   MultiFluidData::CopyCloseNodes copy = MultiFluidData::FALSE;
//   bool printwarning = false;
//   double maxres = -FLT_MAX;
//   int maxreslvl = 1;

//   for (int ilvl = 2; ilvl <= max_level; ++ilvl) {
//     double res = 1.0;
//     double resn = 1.0;
//     double resnm1 = 1.0;

//     int it = 0;
//     while (res > iod.eqs.tc.tm.d2wall.eps && it < iod.eqs.tc.tm.d2wall.maxIts)
//     {
//       resnm1 = resn;
//       dom.computeDistanceLevelNodes(1, tag, ilvl, X, d2wall, resn, dummyPhi, copy);
//       dom.getCommunicator()->globalMax(1, &resn);
//       res = fabs((resn - resnm1) / (resn + resnm1));
//       it++;

//       // dom.getCommunicator()->fprintf(stderr, "Residual = %e @ iteration %d\n", res, it-1);
//     }

//     if (res > iod.eqs.tc.tm.d2wall.eps)
//     {
//       printwarning = true;
//       if (res > maxres)
//       {
//         maxres = res;
//         maxreslvl = ilvl;
//       }
//     }

//     // dom.getCommunicator()->fprintf(stderr, "Wall distance performed maximum %d iterations at level %d of %d\n\n", --it, ilvl, max_level);
//   }

//   if (printwarning)
//     dom.getCommunicator()->fprintf(stderr,
//                                     "*** Warning: Iterative distance to wall computation (Max residual: %e at level: %d, target: %e)\n",
//                                     maxres, maxreslvl, iod.eqs.tc.tm.d2wall.eps);
// }

//------------------------------------------------------------------------------

// predictor function to check if wall distance update is required
// error code key:  0 = no tolerances exceeded, do not update d2wall
//                  1 = change in d = 0, full update d
//                  2 = delta tolerance exceeded, full update d
//                  3 = not enough timestep info to initialize predictors, full update d
template <int dimLS, int dim>
int ReinitializeDistanceToWall<dimLS,dim>::UpdatePredictorsCheckTol(
  DistLevelSetStructure &LSS, DistGeoState &distGeoState, const double t)
{
  int update = 0;

  if (predictorTime[1] < 0)
    update = 3;
  else if (predictorTime[2] < 0)
    update = 2;

  else { // compute distance predictions
    double dtnm2, dtnm1, dtn;
    dtnm2 = predictorTime[1]-predictorTime[2];
    dtnm1 = predictorTime[0]-predictorTime[1];
    dtn = t-predictorTime[0];

    double dtnodtnm1, dtnodtnm2, dtnm1odtnm2, dtnpdtnm1odtnm1pdtnm2;

    double d2wp, delta, meandelta = 0.0;  // RMS error
    double dd, maxdd = -1.0;              // min error

#ifdef PREDICTOR_DEBUG
    dom.getCommunicator()->fprintf(stderr,"\nUpdating predictors...\n");
    double maxdelta = -1.0, maxd2w = -1.0, mind2w = 1.0e10, meandd2w = 0.0;
    double maxdeltaloc;
    // int nvi = 0;
#endif

#pragma omp parallel for
    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
      for (int k = 0; k < d2wall(iSub).len; k++) {
        if (d2wall(iSub)[k][0] > 0.0) {
          if (!(d2wnm1(iSub)[k][0] > 0.0 && d2wnm2(iSub)[k][0] > 0.0)) continue;

          // predictor distance
          dtnodtnm1 = dtn/dtnm1;
          dtnodtnm2 = dtn/dtnm2;
          dtnm1odtnm2 = dtnm1/dtnm2;
          dtnpdtnm1odtnm1pdtnm2 = (dtn+dtnm1)/(dtnm1+dtnm2);
          d2wp = (1.0+dtnodtnm1*(1.0+dtnpdtnm1odtnm1pdtnm2))*d2wall(iSub)[k][0]
            - dtnodtnm1*(dtnodtnm2+dtnm1odtnm2+1.0)*d2wnm1(iSub)[k][0]
            + dtnodtnm2*dtnpdtnm1odtnm1pdtnm2*d2wnm2(iSub)[k][0];
          // d2wp = (1.0+dtnodtnm1*(dtn+2.0*dtnm1+dtnm2)/(dtnm1+dtnm2))*d2wall(iSub)[k][0]
          // -(dtn/dtnm1)*(dtn+dtnm1+dtnm2)/dtnm2*d2wnm1(iSub)[k][0]
          // +(dtn/dtnm2)*(dtn+dtnm1)/(dtnm1+dtnm2)*d2wnm2(iSub)[k][0];

          // // 1/d^2 delta
          // delta = 2.0*fabs(d2wp-d2wall(iSub)[k][0])/d2wall(iSub)[k][0];

          // new delta definition based on SA source error
          delta = fabs(SAsensi(iSub)[k]*(d2wp-d2wall(iSub)[k][0]));
          meandelta += delta*delta;
          dd = fabs(d2wp-d2wall(iSub)[k][0])/d2wall(iSub)[k][0];
          maxdd  = max(dd,maxdd);

#ifdef PREDICTOR_DEBUG
          // for viewing maximum predicted change and other values
          if (delta>maxdelta) {
            maxdelta = delta;
            maxdeltaloc = d2wall(iSub)[k][0];
          }
          mind2w = min(mind2w,d2wall(iSub)[k][0]);
          maxd2w = max(maxd2w,d2wall(iSub)[k][0]);
          meandd2w += dd*dd;
#endif
        }
        else if (LSS(iSub).isActive(0.0,k)) { // ghost real transition

#ifdef PREDICTOR_DEBUG
          fprintf(stderr,"\nGhost to real transition: new d2w = %e (old d2w = %e)\n",
            LSS(iSub).distToInterface(0.0,k),d2wall(iSub)[k][0]);
          if (!(LSS(iSub).isNearInterface(0.0,k))) fprintf(stderr,"Problem: updated node has no exact distance!\n\n");
#endif
          d2wall(iSub)[k][0] = LSS(iSub).distToInterface(0.0,k);
          distGeoState(iSub).getDistanceToWall()[k] = d2wall(iSub)[k][0];
        }
      }
    }

    // RMS error
    dom.getCommunicator()->globalSum(1,&meandelta);
    meandelta = sqrt(meandelta/nPredTot);

#ifdef PREDICTOR_DEBUG
    // max and mean error outputs to view (for testing)
    dom.getCommunicator()->globalMax(1,&maxd2w);
    dom.getCommunicator()->globalMin(1,&mind2w);
    dom.getCommunicator()->fprintf(stderr,"Max predictor distance = %e, min = %e\n",
      maxd2w, mind2w);

    // dom.getCommunicator()->globalMax(1,&maxdelta);
    // dom.getCommunicator()->fprintf(stderr,"Max error = %e (vs. tol = %e)\n",
    //   maxdelta, tolmax);

    struct {
        double errval;
        int d2wval;
    } in, out;
    in.errval = maxdelta;
    in.d2wval = (int) (maxdeltaloc*1.0e7);
    MPI_Comm comm = dom.getCommunicator()->getMPIComm();
    MPI_Allreduce(&in, &out, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, comm);
    dom.getCommunicator()->fprintf(stderr,"Max error = %e at d2w = %e\n",
      out.errval, ((double) out.d2wval)*1.0e-7);

    dom.getCommunicator()->fprintf(stderr,"RMS error = %e  (vs. tol = %e)\n", meandelta, tolmax);

    dom.getCommunicator()->globalSum(1,&meandd2w);
    meandd2w = sqrt(meandd2w/nPredTot);
    dom.getCommunicator()->fprintf(stderr,"RMS delta d2w = %e\n", meandd2w);

    dom.getCommunicator()->fprintf(stderr,"\n");

    // dom.getCommunicator()->globalSum(1,&nvi);
    // dom.getCommunicator()->fprintf(stderr,"Total %d violating predictors\n",nvi);
#endif

    if (meandelta > tolmax)
      update = 2;
    else {
      dom.getCommunicator()->globalMax(1,&maxdd);
      if (maxdd < tolmin)
        update = 1;
    }
  }

  // dom.getCommunicator()->globalMax(1,&update);
  return update;
}

//------------------------------------------------------------------------------

// create new predictors and populate
template <int dimLS, int dim>
void ReinitializeDistanceToWall<dimLS,dim>::ReinitializePredictors(int update,
                                                               const double t,
                                                               DistLevelSetStructure *LSS,
                                                               DistSVec<double, 3> &X)
{
  // pointers to relevant data
  FemEquationTerm *fet = spaceOp.getFemEquationTerm();
  DistSVec<double,dim>& V = *spaceOp.getCurrentPrimitiveVector();
  DistNodalGrad<dim, double>& ngrad = *spaceOp.getDistNodalGrad(V);
  DistSVec<double,dim>& dVdx = ngrad.getX();
  DistSVec<double,dim>& dVdy = ngrad.getY();
  DistSVec<double,dim>& dVdz = ngrad.getZ();

  // non-loop dependent quantities
  const double sixth = 1.0/6.0;
  // double cv1 = iod.eqs.tc.tm.sa.cv1;
  // double cv13 = cv1*cv1*cv1;
  // double cb1 = iod.eqs.tc.tm.sa.cb1;
  // double cb2 = iod.eqs.tc.tm.sa.cb2;
  // double cw2 = iod.eqs.tc.tm.sa.cw2;
  // double cw3 = iod.eqs.tc.tm.sa.cw3;
  // double cw36 = cw3*cw3*cw3*cw3*cw3*cw3;
  // double oosigma = 1.0 / iod.eqs.tc.tm.sa.sigma;
  // double ooK2 = 1.0 / (iod.eqs.tc.tm.sa.vkcst*iod.eqs.tc.tm.sa.vkcst);
  // double cw1 = cb1 * ooK2 + (1.0+cb2) * oosigma;
  // double rlim = 2.0;

  // double ooreynolds_mu = sa->get_ooreynolds_mu();
  // cw1 *= ooreynolds_mu;
  // oosigma *= ooreynolds_mu;

  // instead, pull from sa or des object (as applicable)
  double cb1, cb2, cw1, cw2, cw36, cv13, oosigma, ooK2, rlim, opcw36sixth, ooreynolds_mu;
  ViscoFcn *viscoFcn;
  VarFcn *varFcn;
  FemEquationTermSA* sa = dynamic_cast<FemEquationTermSA*>(fet);
  FemEquationTermDES* des = dynamic_cast<FemEquationTermDES*>(fet);
  if (sa) {
    cb1 = sa->cb1;
    cb2 = sa->cb2;
    cw1 = sa->cw1;
    cw2 = sa->cw2;
    cw36 = sa->cw3_pow6;
    cv13 = sa->cv1_pow3;
    oosigma = sa->oosigma;
    ooK2 = sa->oovkcst2;
    rlim = sa->rlim;
    opcw36sixth = sa->opcw3_pow;
    ooreynolds_mu = sa->oorey;

    varFcn = sa->getVarFcn();
    viscoFcn = sa->getViscoFcn();
  }
  else if (des) {
    cb1 = des->cb1;
    cb2 = des->cb2;
    cw1 = des->cw1;
    cw2 = des->cw2;
    cw36 = des->cw3_pow6;
    cv13 = des->cv1_pow3;
    oosigma = des->oosigma;
    ooK2 = des->oovkcst2;
    rlim = des->rlim;
    opcw36sixth = des->opcw3_pow;
    ooreynolds_mu = des->oorey;

    varFcn = des->getVarFcn();
    viscoFcn = des->getViscoFcn();
  }
  else {
    fprintf(stderr,"*** Error: Cannot construct a NavierStokesTerm in ReinitializeDistanceToWall\n");
    exit(-1);
  }

  // allocate loop dependent quantities
  double rho, T, mul, nul, nutilde, nutilde2;
  double chi, chi3, fv1, fv2, ood2w2, zz;
  double dnutildedxj[3], drhodxj[3];
  double om12, om23, om31, Omega, Stilde;
  double rr, rr5, gg, gg6, oog6cw36, cw36g6sixth, fw;
  double AA, BB, CC, DD, SAterm, drdStilde, dgdr, dfwdg, dStildedd, Lambda, dSAterm;

  nPredTot = 0;

#ifdef PREDICTOR_DEBUG
  double meandSA = 0.0;
#endif

#pragma omp parallel for
  for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
    for (int k = 0; k < d2wall(iSub).len; k++) {
      if (!(d2wall(iSub)[k][0]>0.0 && d2wnm1(iSub)[k][0]>0.0 && d2wnm2(iSub)[k][0]>0.0)) continue;

      // source term and sensitivity calculation
      rho = V(iSub)[k][0];
      T = varFcn->computeTemperature(V(iSub)[k]);
      mul = viscoFcn->compute_mu(T);
      nul = mul/rho;
      nutilde = max(V(iSub)[k][5],0.0);
      nutilde2 = nutilde*nutilde;

      chi = nutilde/nul;
      chi3 = chi*chi*chi;
      fv1 = chi3/(chi3+cv13);
      fv2 = 1.0-chi/(1.0+chi*fv1);

      ood2w2 = 1.0/(d2wall(iSub)[k][0]*d2wall(iSub)[k][0]);
      zz = ooreynolds_mu*ooK2*nutilde*ood2w2;

      // nodal gradients
      dnutildedxj[0] = dVdx(iSub)[k][5];
      dnutildedxj[1] = dVdy(iSub)[k][5];
      dnutildedxj[2] = dVdz(iSub)[k][5];
      drhodxj[0] = dVdx(iSub)[k][0];
      drhodxj[1] = dVdy(iSub)[k][0];
      drhodxj[2] = dVdz(iSub)[k][0];

      // vorticity
      om12 = dVdy(iSub)[k][1] - dVdx(iSub)[k][2];
      om23 = dVdz(iSub)[k][2] - dVdy(iSub)[k][3];
      om31 = dVdx(iSub)[k][3] - dVdz(iSub)[k][1];
      Omega = sqrt(om12*om12 + om23*om23 + om31*om31);
      Stilde = max(Omega+zz*fv2,1.0e-12);

      rr = min(zz/Stilde,rlim);
      rr5 = rr*rr*rr*rr*rr;
      gg = rr+cw2*rr*(rr5-1.0);
      gg6 = gg*gg*gg*gg*gg*gg;
      oog6cw36 = 1.0/(gg6+cw36);
      cw36g6sixth = opcw36sixth*pow(oog6cw36,sixth);
      fw = gg*cw36g6sixth;

      // source term
      AA = rho*cb2*oosigma*(dnutildedxj[0]*dnutildedxj[0]
        + dnutildedxj[1]*dnutildedxj[1] + dnutildedxj[2]*dnutildedxj[2]);
      BB = rho*nutilde*cb1*Stilde;
      CC = -rho*cw1*fw*nutilde2*ood2w2;
      DD = -(nul+nutilde)*oosigma*(dnutildedxj[0]*drhodxj[0]
        + dnutildedxj[1]*drhodxj[1] + dnutildedxj[2]*drhodxj[2]);
      SAterm = AA + BB + CC + DD;

      // // alternatively, only look at sensitivity of nutilde production and destruction
      // BB = rho*nutilde*cb1*Stilde;
      // CC = -rho*cw1*fw*nutilde2*ood2w2;
      // SAterm = BB + CC;

      // sensitivity
      drdStilde = (rr>rlim)?-zz/(Stilde*Stilde):0.0;
      dgdr = 1.0+cw2*(6.0*rr5-1.0);
      dfwdg = cw36*oog6cw36*cw36g6sixth;
      dStildedd = -2.0*fv2*zz/d2wall(iSub)[k][0];
      Lambda = cb1*nutilde-cw1*nutilde2*ood2w2*dfwdg*dgdr*drdStilde;
      dSAterm = rho*(2.0*cw1*fw*nutilde2*(ood2w2*d2wall(iSub)[k][0]) + Lambda*dStildedd);

#ifdef PREDICTOR_DEBUG
      if (SAsensi(iSub)[k]>0.0)
        meandSA += ((dSAterm/SAterm-SAsensi(iSub)[k])/SAsensi(iSub)[k])
          *((dSAterm/SAterm-SAsensi(iSub)[k])/SAsensi(iSub)[k]);
#endif

      // store sensitivities and predictor count
      SAsensi(iSub)[k] = dSAterm/SAterm;
      nPredTot++;
    }
  }
  dom.getCommunicator()->globalSum(1,&nPredTot);

#ifdef PREDICTOR_DEBUG
  dom.getCommunicator()->globalSum(1,&meandSA);
  meandSA = sqrt(meandSA/nPredTot);
  dom.getCommunicator()->fprintf(stderr,"RMS delta SA sensitivity = %e. ", meandSA);
  dom.getCommunicator()->fprintf(stderr, "Initialized %d predictors\n\n", nPredTot);
#endif
}

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

  // // compare to iterative
  // int iterativeTemp = iod.eqs.tc.tm.d2wall.maxIts;
  // iod.eqs.tc.tm.d2wall.maxIts = 100;
  // dom.getCommunicator()->fprintf(stderr,
  //   "Comparing specified wall distance computation to iterative method with 100 maximum iterations.\n");
  // DistanceToClosestPointOnMovingStructure(LSS, X, distGeoState);
  // GetLevelsFromInterfaceAndMarchForward(LSS, X, distGeoState);
  // iod.eqs.tc.tm.d2wall.maxIts = iterativeTemp;

  // compare to noniterative
  dom.getCommunicator()->fprintf(stderr,
    "Comparing specified wall distance computation to noniterative method.\n");
  PseudoFastMarchingMethod(LSS,X,0);
  // PseudoFastMarchingMethodFEM(LSS,X);

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

// sjg, 02/2017: compute wall distance change from previous call for testing
template <int dimLS, int dim>
void ReinitializeDistanceToWall<dimLS,dim>::ComputePercentChange(DistLevelSetStructure &LSS, DistGeoState &distGeoState)
{
  int nSub = dom.getNumLocSub();
  double **d2wChange;
  int nDofs[nSub];
  d2wChange = new double *[nSub];
  double localDelta = 0.0;
  DistSVec<double, 1> d2wallNew = d2wall;
  DistSVec<double, 1> d2wallPrev(dom.getNodeDistInfo());

  // Check if first timestep and return if so (nothing to compare to)
  int iSub, i;
#pragma omp parallel for
  for (iSub = 0; iSub < nSub; ++iSub)
  {
    for (i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i)
    {
      if (distGeoState(iSub).getDistanceToWall()[i] > 0.0)
      {
        break;
      }
    }
    if (i < distGeoState(iSub).getDistanceToWall().size())
    {
      break;
    }
  }
  if (iSub == nSub && i == distGeoState(iSub - 1).getDistanceToWall().size())
  {
    // dom.getCommunicator()->fprintf(stderr,"First time step, nothing to compare to for percent change.\n");
    return;
  }

// Populate for storing old values
#pragma omp parallel for
  for (int iSub = 0; iSub < nSub; ++iSub)
    for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i)
      d2wallPrev(iSub)[i][0] = distGeoState(iSub).getDistanceToWall()[i];

#pragma omp parallel for
  for (int iSub = 0; iSub < nSub; ++iSub)
  {
    d2wChange[iSub] = new double[4];
    d2wChange[iSub][0] = 0.0;
    d2wChange[iSub][1] = 0.0;
    d2wChange[iSub][2] = 0.0;
    d2wChange[iSub][3] = 0.0;

    nDofs[iSub] = 0;
    for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i)
    {
      if (LSS(iSub).isActive(0.0, i) && d2wallPrev(iSub)[i][0] > 0.0)
      {
        nDofs[iSub]++;
        localDelta = fabs(d2wallNew(iSub)[i][0] - d2wallPrev(iSub)[i][0]);
        d2wChange[iSub][0] += localDelta / d2wallPrev(iSub)[i][0];
        d2wChange[iSub][1] += localDelta;
        if (localDelta > d2wChange[iSub][2])
        {
          d2wChange[iSub][2] = localDelta;
          d2wChange[iSub][3] = d2wallPrev(iSub)[i][0];
        }
      }
    }
  }
  for (int iSub = 1; iSub < nSub; ++iSub)
  {
    nDofs[0] += nDofs[iSub];
    d2wChange[0][0] += d2wChange[iSub][0];
    d2wChange[0][1] += d2wChange[iSub][1];
    if (d2wChange[iSub][2] > d2wChange[0][2])
    {
      d2wChange[0][2] = d2wChange[iSub][2];
      d2wChange[0][3] = d2wChange[iSub][3];
    }
  }

  // Communicate across all processes to find global sums/max
  dom.getCommunicator()->globalSum(1, nDofs);
  dom.getCommunicator()->globalSum(2, d2wChange[0]);
  MPI_Comm comm = dom.getCommunicator()->getMPIComm();
  MPI_Allreduce(&d2wChange[0][2], &d2wChange[0][2], 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, comm);

  d2wChange[0][0] /= nDofs[0];
  d2wChange[0][1] /= nDofs[0];

  dom.getCommunicator()->fprintf(stderr, "Average relative and absolute d2wall change since last call: %12.8e, %12.8e\n", d2wChange[0][0], d2wChange[0][1]);
  dom.getCommunicator()->fprintf(stderr, "Maximum absolute d2wall change since last call: %12.8e at %12.8e\n\n", d2wChange[0][2], d2wChange[0][3]);

  for (int iSub = 0; iSub < nSub; iSub++)
    delete[] d2wChange[iSub];
  delete[] d2wChange;
  return;
}

//------------------------------------------------------------------------------

// sjg, 04/2017: print wall distance at intersected edges
template <int dimLS, int dim>
void ReinitializeDistanceToWall<dimLS,dim>::PrintIntersectedValues(DistLevelSetStructure *LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState)
{
  dom.getCommunicator()->barrier();
  dom.getCommunicator()->fprintf(stderr, "Printing distance at intersected edges\n\n");

  int nSub = dom.getNumLocSub();
  double **vtag;
  vtag = new double*[nSub];

  for (int iSub = 0; iSub < nSub; ++iSub) {
    LevelSetStructure *locLSS = &((*LSS)(iSub));
    vtag[iSub] = new double[distGeoState(iSub).getDistanceToWall().size()];
    std::memset(vtag[iSub], 0, sizeof(double)*distGeoState(iSub).getDistanceToWall().size());

    int(*ptrEdge)[2] = (*dom.getSubDomain()[iSub]).getEdges().getPtr();
    for (int l = 0; l < (*dom.getSubDomain()[iSub]).getEdges().size(); ++l)
    {
      if (locLSS->edgeIntersectsStructure(0, l))
      {
        int i = ptrEdge[l][0];
        int j = ptrEdge[l][1];
        bool iActive = locLSS->isActive(0.0,i);
        bool jActive = locLSS->isActive(0.0,j);

        if (iActive && !vtag[iSub][i]) {
          Vec3D pt = X(iSub)[i];
          vtag[iSub][i] = 1;
          fprintf(stderr, "d2w(cpu %3d, lsubdom %3d, node %5d at %12.8e,%12.8e,%12.8e) = %12.8e\n",
            dom.getCommunicator()->cpuNum(),iSub,i,pt[0],pt[1],pt[2],d2wall(iSub)[i][0]);
        }
        if (jActive && !vtag[iSub][i]) {
          Vec3D pt = X(iSub)[j];
          vtag[iSub][j] = 1;
          fprintf(stderr, "d2w(cpu %3d, lsubdom %3d, node %5d at %12.8e,%12.8e,%12.8e) = %12.8e\n",
            dom.getCommunicator()->cpuNum(),iSub,j,pt[0],pt[1],pt[2],d2wall(iSub)[j][0]);
        }
      }
    }
    delete[] vtag[iSub];
  }
  delete[] vtag;

  dom.getCommunicator()->barrier();
  dom.getCommunicator()->fprintf(stderr, "\n");
}

//------------------------------------------------------------------------------

template class ReinitializeDistanceToWall<1,5>;
template class ReinitializeDistanceToWall<1,6>;
template class ReinitializeDistanceToWall<1,7>;
template class ReinitializeDistanceToWall<2,5>;
template class ReinitializeDistanceToWall<2,6>;
template class ReinitializeDistanceToWall<2,7>;