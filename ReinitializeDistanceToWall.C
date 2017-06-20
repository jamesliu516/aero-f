#include <ReinitializeDistanceToWall.h>
#include <LevelSet/LevelSetStructure.h>
#include <Domain.h>
#include <DistVector.h>

// sjg, 02/2017: flags for testing and error calculations, and wall distance predictors
  // #define DELTA_CHECK
  // #define ERROR_CHECK
  // #define PRINT_INTERSECT

  #define USE_PREDICTORS
  #define PREDICTOR_DEBUG

//------------------------------------------------------------------------------

template <int dimLS>
ReinitializeDistanceToWall<dimLS>::ReinitializeDistanceToWall(IoData &ioData, Domain &domain)
    : iod(ioData), dom(domain), d2wall(domain.getNodeDistInfo()), tag(domain.getNodeDistInfo()),
    sortedNodes(domain.getNodeDistInfo()), countreinits(0), nPredTot(0),
    d2wnm1(domain.getNodeDistInfo()), d2wnm2(domain.getNodeDistInfo())
{
  int nSub = dom.getNumLocSub();
  nSortedNodes = new int[nSub];
  nActiveNodes = new int[nSub];
  firstCheckedNode = new int[nSub];

  predictorTime = new double[3];
  predictorTime[0] = -1.0;
  predictorTime[1] = -1.0;
  predictorTime[2] = -1.0;

  // predictor tolerances, to be read from parser eventually
  tolmax = 1.0e-2, tolmin = 1.0e-9;
}

//------------------------------------------------------------------------------

template <int dimLS>
ReinitializeDistanceToWall<dimLS>::~ReinitializeDistanceToWall()
{

#ifdef PREDICTOR_DEBUG
  dom.getCommunicator()->fprintf(stderr,"\nWall distance computer: total full domain reinitializations = %d\n\n",countreinits);
#endif

  delete[] nSortedNodes;
  delete[] nActiveNodes;
  delete[] firstCheckedNode;
  delete[] predictorTime;
}

//------------------------------------------------------------------------------

template <int dimLS>
void ReinitializeDistanceToWall<dimLS>::ComputeWallFunction(DistLevelSetStructure &LSS,
                                                            DistSVec<double, 3> &X,
                                                            DistGeoState &distGeoState,
                                                            double t)
{

#ifdef USE_PREDICTORS
  // update predictors and check tolerances, receiving global maximum
  int update = UpdatePredictorsCheckTol(t, LSS, distGeoState);
#else
  int update = 3;
#endif

  // predicted distance change is zero, so correct predictor to exact values
  if (update == 1) {
    // store with previous exact updates
    predictorTime[2] = predictorTime[1];
    predictorTime[1] = predictorTime[0];
    predictorTime[0] = t;
    d2wnm2 = d2wnm1;
    d2wnm1 = d2wall;

    ReinitializePredictors(update, t, &LSS, X);
  }
  else if (update > 1) {
    // compute distance in full domain
#ifdef PREDICTOR_DEBUG
    countreinits++;
    dom.getCommunicator()->fprintf(stderr, "Performing a full domain wall distance update (update = %d)!\n",update);
#endif

    // store with previous exact updates
    predictorTime[2] = predictorTime[1];
    predictorTime[1] = predictorTime[0];
    predictorTime[0] = t;
    d2wnm2 = d2wnm1;
    d2wnm1 = d2wall;

    // update d2wall exactly
    if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::ITERATIVE) {
      // InitializeWallFunction(LSS, X, distGeoState, t);
      // GetLevelsFromInterfaceAndMarchForward(LSS, X, distGeoState);

      // new iterative implementation
      IterativeMethodUpdate(LSS, X);
    }
    else if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::NONITERATIVE) {
      PseudoFastMarchingMethod(LSS, X, 0);
    }
    // else if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::HYBRID) {
    //   int iterativeLevel = 0;

    //   if (iod.eqs.tc.tm.d2wall.iterativelvl > 1) {
    //     // InitializeWallFunction(LSS, X, distGe USE_PREDICTORS
    //     // GetLevelsFromInterfaceAndMarchForward(LSS, X, distGeoState);

    //     IterativeMethodUpdate(LSS, X);
    //     iterativeLevel = iod.eqs.tc.tm.d2wall.iterativelvl;
    //   }
    //   PseudoFastMarchingMethod(LSS, X, iterativeLevel);
    // }
    else {
      fprintf(stderr, " *** Error ***, Unknown wall distance method\n");
      exit(1);
    }

    // reinitialize predictors
    if (update == 2)
      ReinitializePredictors(update, t, &LSS, X);

#ifdef DELTA_CHECK
    ComputePercentChange(LSS, distGeoState); // sjg, 02/2017: percent change since last call for testing
#endif

#pragma omp parallel for
    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
      for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i) {
        distGeoState(iSub).getDistanceToWall()[i] = d2wall(iSub)[i][0];
      }
    }

#ifdef PRINT_INTERSECT
    PrintIntersectedValues(&LSS, X, distGeoState); // sjg, 04/2017: wall distance change in time at intersected edges
#endif
  }
#ifdef PREDICTOR_DEBUG
  else {
    dom.getCommunicator()->fprintf(stderr, "SUCCESS: Skipped wall distance calculation (update = %d)!\n\n",update);
  }
  // dom.getCommunicator()->fprintf(stderr, "Times: tn = %e, tnm1 = %e, tnm2 = %e\n\n",predictorTime[0],predictorTime[1],predictorTime[2]);
#endif

#ifdef ERROR_CHECK
  ComputeExactErrors(LSS, X, distGeoState);
#endif

  return;
}

//------------------------------------------------------------------------------

// sjg, 04/2017: predictor function to check if wall distance update is required
// error code key:  0 = no tolerances exceeded, do not update d2wall
//                  1 = change in d = 0, reinitialize d @ current predictors
//                  2 = delta tolerance exceeded, reinitialize and update d in flow domain
//                  3 = not enough timestep info to initialize predictors
template <int dimLS>
int ReinitializeDistanceToWall<dimLS>::UpdatePredictorsCheckTol(double t,
  DistLevelSetStructure &LSS, DistGeoState &distGeoState)
{
  int update = 0;

  if (predictorTime[1] < 0)
    update = 3;
  else if (predictorTime[2] < 0)
    update = 2;
  else {
    // compute distance predictions
#ifdef PREDICTOR_DEBUG
    dom.getCommunicator()->fprintf(stderr,"\nUpdating predictors...\n");
#endif

    double dtnm2, dtnm1, dtn;
    dtnm2 = predictorTime[1]-predictorTime[2];
    dtnm1 = predictorTime[0]-predictorTime[1];
    dtn = t-predictorTime[0];

    double maxdd = -1.0;
    double delta, d2wp, deltatmp, deltam, meandel = 0.0;
    int k;

#ifdef PREDICTOR_DEBUG
    double maxdel = -1.0, maxd2w = -1.0, mind2w = 1.0e16;
    // int nvi = 0;
#endif

#pragma omp parallel for
    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
      for (int k = 0; k < tag(iSub).len; k++) {
        if (d2wall(iSub)[k][0] > 0.0) {
          if (d2wnm1(iSub)[k][0] > 0.0 && d2wnm2(iSub)[k][0] > 0.0) {

            // d2wp = (1.0+dtn/dtnm1*(dtn+2.0*dtnm1+dtnm2)/(dtnm1+dtnm2))*d2wnm1(iSub)[k][0]
            // -(dtn/dtnm1)*(dtn+dtnm1+dtnm2)/dtnm2*d2wnm2(iSub)[k][0]
            // +(dtn/dtnm2)*(dtn+dtnm1)/(dtnm1+dtnm2)*d2wnm3(iSub)[k][0];
            d2wp = (1.0+dtn/dtnm1*(dtn+2.0*dtnm1+dtnm2)/(dtnm1+dtnm2))*d2wall(iSub)[k][0]
            -(dtn/dtnm1)*(dtn+dtnm1+dtnm2)/dtnm2*d2wnm1(iSub)[k][0]
            +(dtn/dtnm2)*(dtn+dtnm1)/(dtnm1+dtnm2)*d2wnm2(iSub)[k][0];

            delta = 2.0*fabs(d2wp-d2wall(iSub)[k][0])/d2wall(iSub)[k][0];
            meandel += delta*delta;
#ifdef PREDICTOR_DEBUG
            // for viewing maximum predicted change and other values
            maxdel = max(delta,maxdel);
            mind2w = min(mind2w,d2wall(iSub)[k][0]);
            maxd2w = max(maxd2w,d2wall(iSub)[k][0]);
#endif

      //       if (delta > tolmax) {
      //         update = 2;
      // // #ifdef PREDICTOR_DEBUG
      //         fprintf(stderr,
      //           "Tolerance exceeded (2delta/d = %e > %e @ cpu %d, node %d\n    d2wp = %e, d2wall(tnm1) = %e\n    d2wnm1 = %e, d2wnm2 = %e)\n",
      //           delta,tolmax,dom.getCommunicator()->cpuNum(),j,d2wp,d2wall(iSub)[k][0],
      //           d2wnm1(iSub)[k][0],d2wnm2(iSub)[k][0]);
      //         nvi++;
      // // #else
      //         goto exitlbl;
      // // #endif
      //       }
            // else {
              deltam = fabs(d2wp-d2wall(iSub)[k][0])/d2wall(iSub)[k][0];
              maxdd  = max(deltam,maxdd);
            // }
          }
        }
        else if (LSS(iSub).isActive(0.0,k)) {
          // ghost real transition

#ifdef PREDICTOR_DEBUG
          // debugging!!
          fprintf(stderr,"\nGhost to real transition: new d2w = %e (old d2w = %e)\n\n",
            LSS(iSub).distToInterface(0.0,k),d2wall(iSub)[k][0]);
#endif
          d2wall(iSub)[k][0] = LSS(iSub).distToInterface(0.0,k);
          distGeoState(iSub).getDistanceToWall()[k] = d2wall(iSub)[k][0];
        }
      }
    }

    // exitlbl:
    // dom.getCommunicator()->globalMax(1,&update);

    // RMS error
    dom.getCommunicator()->globalSum(1,&meandel);
    meandel = sqrt(meandel/nPredTot);

#ifdef PREDICTOR_DEBUG
    // max and mean error outputs to view (for testing)
    dom.getCommunicator()->globalMax(1,&maxdel);
    dom.getCommunicator()->fprintf(stderr,"Max (2deltad/d) = %e (vs. tol = %e)\n",
      maxdel, tolmax);

    dom.getCommunicator()->fprintf(stderr,"RMS (2deltad/d) = %e\n", meandel);

    dom.getCommunicator()->globalMax(1,&maxd2w);
    dom.getCommunicator()->globalMin(1,&mind2w);
    dom.getCommunicator()->fprintf(stderr,"Max predictor distance = %e, min = %e\n\n",
      maxd2w, mind2w);
#endif

    if (meandel > tolmax)
      update = 2;
    else {
    // if (update < 2) {
      dom.getCommunicator()->globalMax(1,&maxdd);
      if (maxdd < tolmin)
        update = 1;
    }

  // #ifdef PREDICTOR_DEBUG
  //     dom.getCommunicator()->globalSum(1,&nvi);
  //     dom.getCommunicator()->fprintf(stderr,"Total %d violating predictors\n",nvi);
  // #endif
  }

  // dom.getCommunicator()->globalMax(1,&update);
  return update;
}

//------------------------------------------------------------------------------

// sjg, 04/2017: reinitialize predictors
template <int dimLS>
void ReinitializeDistanceToWall<dimLS>::ReinitializePredictors(int update,
                                                               double t,
                                                               DistLevelSetStructure *LSS,
                                                               DistSVec<double, 3> &X)
{
  int nSub = dom.getNumLocSub();
  // update distance at current predictors exactly
  if (update == 1) {

#ifdef PREDICTOR_DEBUG
    dom.getCommunicator()->fprintf(stderr, "Reinitializing due to flat ...\n");
#endif

    sortedNodes = -1;
    tag = -1;
    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
      nSortedNodes[iSub] = 0;
      nActiveNodes[iSub] = 0;
      firstCheckedNode[iSub] = 0;
    }
    dom.pseudoFastMarchingMethod<1>(tag, X, d2wall, 1, 0,
      sortedNodes, nSortedNodes, nActiveNodes, firstCheckedNode, LSS);
  }

  // create new predictors and populate
  else if (update == 2) {

#ifdef PREDICTOR_DEBUG
    dom.getCommunicator()->fprintf(stderr, "Reinitializing predictors...");
#endif

    nPredTot = 0;
#pragma omp parallel for
    for (int iSub = 0; iSub < nSub; ++iSub) {
      for (int i = 0; i < tag(iSub).len; i++)
        nPredTot += (d2wall(iSub)[i][0]>0.0 && d2wnm1(iSub)[i][0]>0.0 && d2wnm2(iSub)[i][0]>0.0);
    }
    dom.getCommunicator()->globalSum(1,&nPredTot);

#ifdef PREDICTOR_DEBUG
    dom.getCommunicator()->fprintf(stderr, "initialized %d predictors\n\n", nPredTot);
#endif
  }
  else {
      fprintf(stderr, " *** Error ***, Wall distance predictor update case is invalid\n");
      exit(1);
  }
}

//------------------------------------------------------------------------------

// The following is an adaptation of the Fast Marching Method to Embedded Turbulent computation.
// Adam 2012.09
template <int dimLS>
int ReinitializeDistanceToWall<dimLS>::PseudoFastMarchingMethod(
    DistLevelSetStructure &LSS, DistSVec<double, 3> &X, int iterativeLevel)
{
  sortedNodes = -1;
  // if (iterativeLevel == 0)
  // {
    d2wall = 1.0e10;
    tag = -1;
  // }
  int isDone = 0;
  int nSub = dom.getNumLocSub();

  int level = iterativeLevel; // Level 0 (inActive nodes) and 1 (Embedded surface neighbors)
  while (isDone == 0)
  { // Tag and update d at every level
    dom.pseudoFastMarchingMethod<1>(tag, X, d2wall, level, 0, sortedNodes, nSortedNodes, nActiveNodes, firstCheckedNode, &LSS);
    // I don't think it is a good idea to OMP parallelize this loop. nSub should be small, though!
    isDone = 1;
    for (int iSub = 0; iSub < nSub; ++iSub) {
      if (nSortedNodes[iSub] != tag(iSub).len) {
        isDone = 0;
        break;
      }
    }
    dom.getCommunicator()->globalMin(1, &isDone);
    ++level;
  }
  --level;

  // dom.getCommunicator()->fprintf(stderr, "There are %d levels\n", level);
  return level;
}

//------------------------------------------------------------------------------

template <int dimLS>
void ReinitializeDistanceToWall<dimLS>::IterativeMethodUpdate(DistLevelSetStructure &LSS,
                                                              DistSVec<double, 3> &X)
{
  // initial pass (local FMM)
  sortedNodes = -1;
  d2wall = 1.0e10;
  tag = -1;
  int iSub, isDone, level, nSub = dom.getNumLocSub();
  double res0, resnm1, res = 0.0;
  for (iSub = 0; iSub < nSub; ++iSub) {
    level = 0, isDone = 0;
    while (level < 2 && isDone == 0) {
      res += dom.pseudoFastMarchingMethodSerial<1>(iSub, tag, X, d2wall, level, 0,
        sortedNodes, nSortedNodes, nActiveNodes, firstCheckedNode, &LSS);
      isDone = (nSortedNodes[iSub]==tag(iSub).len)?1:0;
      ++level;
    }
    // if (nSortedNodes[iSub] > 0) {
    //   while (isDone == 0) {
    //     dom.pseudoFastMarchingMethodSerial<1>(iSub, tag, X, d2wall, level, 0,
    //       sortedNodes, nSortedNodes, firstCheckedNode, &LSS);
    //     isDone = (nSortedNodes[iSub]==tag(iSub).len)?1:0;
    //     ++level;
    //   }
    // }
  }
  dom.getCommunicator()->globalSum(1, &res);
  res = sqrt(res);

  // // share initialization information across processes
  // dom.pseudoFastMarchingMethodComm<1>(tag, d2wall, sortedNodes, nSortedNodes, nActiveNodes, res);

  /* NOTE: max level and hybrid method no longer possible with full domain
           iteration (enforce instead using a max distance if desired) */
  // int max_level = level-1;
  // if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::HYBRID &&
  //     iod.eqs.tc.tm.d2wall.iterativelvl > 1)
  //   max_level = min(iod.eqs.tc.tm.d2wall.iterativelvl, max_level);

  // communication and iteration loop
  // DistSVec<double, 1> resnm1(dom.getNodeDistInfo());
  // double res[2];
  // res[0] = -1.0; res[1] = -1.0;
  int it = 1, isConverged = 0;

  while (!isConverged) {
    // resnm1 = d2wall;
    resnm1 = res;
    res = 0.0;

    // share initialization information across processes
    dom.pseudoFastMarchingMethodComm<1>(tag, d2wall, sortedNodes, nSortedNodes, nActiveNodes);

    // propagate information outwards from minimum level
    isConverged = 1;
#pragma omp parallel for
    for (iSub = 0; iSub < nSub; ++iSub) {
      if (nSortedNodes[iSub] > 0) {
        level = 1, isDone = 0;
        while (isDone == 0) {
          res += dom.pseudoFastMarchingMethodSerial<1>(iSub, tag, X, d2wall, level, 1,
            sortedNodes, nSortedNodes, nActiveNodes, firstCheckedNode, &LSS);
          // isDone = (nSortedNodes[iSub]==tag(iSub).len)?1:0;
          isDone = (nSortedNodes[iSub]==tag(iSub).len || nSortedNodes[iSub]==0)?1:0;
          ++level;
        }
      }
      else
        isConverged = 0;
    }
    dom.getCommunicator()->globalSum(1, &res);
    res = sqrt(res);

    dom.getCommunicator()->globalMin(1, &isConverged);

    // // share new information across processes
    // dom.pseudoFastMarchingMethodComm<1>(tag, d2wall, sortedNodes, nSortedNodes, nActiveNodes);

    // update residual
    if (isConverged) {
      // res[0] = 0.0;
      // res[1] = 0.0;
      // for (iSub = 0; iSub < nSub; ++iSub) {
      //   for (int i = 0; i < tag(iSub).len; i++) {
      //     if (d2wall(iSub)[i][0] > 0.0) {
      //       // res[0] = max(fabs(resnm1(iSub)[i][0]-d2wall(iSub)[i][0])/resnm1(iSub)[i][0],res[0]);

      //       // 2-norm
      //       // res[0] += (d2wall(iSub)[i][0]-resnm1(iSub)[i][0])*(d2wall(iSub)[i][0]-resnm1(iSub)[i][0]);
      //       // res[1] += resnm1(iSub)[i][0]*resnm1(iSub)[i][0];

      //       // RMS
      //       res[0] += ((d2wall(iSub)[i][0]-resnm1(iSub)[i][0])/resnm1(iSub)[i][0])
      //         *((d2wall(iSub)[i][0]-resnm1(iSub)[i][0])/resnm1(iSub)[i][0]);
      //       res[1]++;
      //     }
      //   }
      // }
      // // dom.getCommunicator()->globalMax(1, res);
      // // res[0] = sqrt(res[0])/sqrt(res[1]);
      // dom.getCommunicator()->globalSum(2, res);
      // // res[0] = sqrt(res[0])/sqrt(res[1]);
      // res[0] = sqrt(res[0]/res[1]);
      // isConverged = (res[0] < iod.eqs.tc.tm.d2wall.eps || it > iod.eqs.tc.tm.d2wall.maxIts);

      res0 = fabs(res-resnm1)/(res+resnm1);
      isConverged = (res0 < iod.eqs.tc.tm.d2wall.eps || it > iod.eqs.tc.tm.d2wall.maxIts);

      // dom.getCommunicator()->fprintf(stderr, "Residual = %e @ iteration %d (isConverged = %d)\n", res[0], it, isConverged);
    }
    it++;
  }

  // dom.getCommunicator()->fprintf(stderr,
  //   "Iterative distance to wall computation: final residual = %e, target = %e @ iteration %d\n",
  //   res[0], iod.eqs.tc.tm.d2wall.eps, it-1);
  dom.getCommunicator()->fprintf(stderr,
    "Iterative distance to wall computation: final residual = %e, target = %e @ iteration %d\n",
    res0, iod.eqs.tc.tm.d2wall.eps, it-1);
}

//------------------------------------------------------------------------------

// template <int dimLS>
// void ReinitializeDistanceToWall<dimLS>::InitializeWallFunction(DistLevelSetStructure &LSS,
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

// template <int dimLS>
// void ReinitializeDistanceToWall<dimLS>::GetLevelsFromInterfaceAndMarchForward(DistLevelSetStructure &LSS,
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
//         min_level = min(min_level, tag(iSub)[i]);
//         if (min_level < 0)
//           goto exitlbl;
//         max_level = max(max_level, tag(iSub)[i]);
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

// sjg, 02/2017: instead of computing error for embedded cylinder, compute relative
// error of iterative and noniterative methods
template <int dimLS>
void ReinitializeDistanceToWall<dimLS>::ComputeExactErrors(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState)
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

  MPI_Comm comm = dom.getCommunicator()->getMPIComm();
  MPI_Allreduce(&errors[0][2], &errors[0][2], 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, comm);
  MPI_Allreduce(&errorsEx[0][2], &errorsEx[0][2], 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, comm);

  errors[0][0] /= nDofs[0];
  errors[0][1] /= nDofs[0];
  errors[0][1] = sqrt(errors[0][1]);
  errorsEx[0][0] /= nDofs[0];
  errorsEx[0][1] /= nDofs[0];
  errorsEx[0][1] = sqrt(errorsEx[0][1]);

  dom.getCommunicator()->fprintf(stderr, "Absolute d2wall Error: %12.8e, %12.8e, %12.8e at %12.8e\n", errors[0][0], errors[0][1], errors[0][2], errors[0][3]);
  dom.getCommunicator()->fprintf(stderr, "Relative d2wall Error: %12.8e, %12.8e, %12.8e at %12.8e\n\n", errorsEx[0][0], errorsEx[0][1], errorsEx[0][2], errorsEx[0][3]);

  for (int iSub = 0; iSub < nSub; iSub++)
    delete[] errors[iSub];
  delete[] errors;
  for (int iSub = 0; iSub < nSub; iSub++)
    delete[] errorsEx[iSub];
  delete[] errorsEx;
  return;
}

//------------------------------------------------------------------------------

template <int dimLS>
void ReinitializeDistanceToWall<dimLS>::PrescribedValues(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState)
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
template <int dimLS>
void ReinitializeDistanceToWall<dimLS>::ComputePercentChange(DistLevelSetStructure &LSS, DistGeoState &distGeoState)
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
  {
    for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i)
      d2wallPrev(iSub)[i][0] = distGeoState(iSub).getDistanceToWall()[i];
  }

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
template <int dimLS>
void ReinitializeDistanceToWall<dimLS>::PrintIntersectedValues(DistLevelSetStructure *LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState)
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
template class ReinitializeDistanceToWall<1>;
template class ReinitializeDistanceToWall<2>;