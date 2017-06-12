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
    : iod(ioData), dom(domain), d2wall(domain.getNodeDistInfo()), tag(domain.getNodeDistInfo()), dummyPhi(domain.getNodeDistInfo()), sortedNodes(domain.getNodeDistInfo()), predictorActive(0),
    d2wallnm1(domain.getNodeDistInfo()), d2wallnm2(domain.getNodeDistInfo()),tnm1(-1.0),tnm2(-1.0),countReinits(0)
{
  int nSub = dom.getNumLocSub();
  nSortedNodes = new int[nSub];
  firstCheckedNode = new int[nSub];

  nPredictors = new int[nSub];
  for (int iSub = 0; iSub<nSub; iSub++)
    nPredictors[iSub] = 0;

  tPredictors = new double[3];
  tPredictors[0] = -1.0;
  tPredictors[1] = -1.0;
  tPredictors[2] = -1.0;

  d2wPredictors = new double **[nSub];
}

//------------------------------------------------------------------------------

template <int dimLS>
ReinitializeDistanceToWall<dimLS>::~ReinitializeDistanceToWall()
{

#ifdef PREDICTOR_DEBUG
  dom.getCommunicator()->fprintf(stderr,"\nWall distance computer: total full domain reinitializations = %d\n\n",countReinits);
#endif

  delete[] nSortedNodes;
  delete[] firstCheckedNode;

  if (predictorActive) {
    for (int iSub = 0; iSub<dom.getNumLocSub(); iSub++) {
      for (int j = 0; j<nPredictors[iSub]; j++) {
        delete[] d2wPredictors[iSub][j];
      }
      delete [] d2wPredictors[iSub];
    }
  }
  delete[] d2wPredictors;

  delete[] tPredictors;
  delete[] nPredictors;
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
  int update = UpdatePredictorsCheckTol(t);
#else
  int update = 3;
#endif

  // predicted distance change is zero, so correct predictor to exact values
  if (update == 1) {
    ReinitializePredictors(update, t, &LSS, X, distGeoState);

    d2wallnm2 = d2wallnm1;
    d2wallnm1 = d2wall;
    tnm2 = tnm1;
    tnm1 = t;
  }

  // compute distance in full domain
  else if (update > 1) {

    countReinits++;
#ifdef PREDICTOR_DEBUG
    dom.getCommunicator()->fprintf(stderr, "Performing a full domain wall distance update (update = %d)!\n",update);
#endif

    // free memory from previous predictors if reinitializing
    if (predictorActive) {
      for (int iSub = 0; iSub<dom.getNumLocSub(); iSub++) {
        for (int j = 0; j<nPredictors[iSub]; j++) {
          delete[] d2wPredictors[iSub][j];
        }
        delete [] d2wPredictors[iSub];
      }
    }

    // update d2wall exactly
    if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::ITERATIVE) {
      // InitializeWallFunction(LSS, X, distGeoState, nPredictors, t);
      // GetLevelsFromInterfaceAndMarchForward(LSS, X, distGeoState);
      int maxlvl = PseudoFastMarchingMethod(LSS, X, distGeoState, 0, nPredictors);
      IterativeMethodUpdate(maxlvl, LSS, X, distGeoState);
    }
    else if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::NONITERATIVE) {
      PseudoFastMarchingMethod(LSS, X, distGeoState, 0, nPredictors);
    }
    else if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::HYBRID) {
      int iterativeLevel = 0;

      if (iod.eqs.tc.tm.d2wall.iterativelvl > 1) {
        // InitializeWallFunction(LSS, X, distGeoState, nPredictors, t);
        // GetLevelsFromInterfaceAndMarchForward(LSS, X, distGeoState);
        int maxlvl = PseudoFastMarchingMethod(LSS, X, distGeoState, 0, nPredictors);
        IterativeMethodUpdate(maxlvl, LSS, X, distGeoState);
        iterativeLevel = iod.eqs.tc.tm.d2wall.iterativelvl;
      }
      PseudoFastMarchingMethod(LSS, X, distGeoState, iterativeLevel, nPredictors);
    }
    else {
      fprintf(stderr, " *** Error ***, Unknown wall distance method\n");
      exit(1);
    }

    // reinitialize predictors
    if (update == 2)
      ReinitializePredictors(update, t, &LSS, X, distGeoState);

    // store previous exact updates
    d2wallnm2 = d2wallnm1;
    d2wallnm1 = d2wall;
    tnm2 = tnm1;
    tnm1 = t;

// sjg, 02/2017: percent change since last call for testing
#ifdef DELTA_CHECK
    ComputePercentChange(LSS, distGeoState);
#endif

#pragma omp parallel for
    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
      for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i) {
        distGeoState(iSub).getDistanceToWall()[i] = d2wall(iSub)[i][0];
      }
    }

// sjg, 04/2017: wall distance change in time at intersected edges
#ifdef PRINT_INTERSECT
    PrintIntersectedValues(&LSS, X, distGeoState);
#endif
  }

#ifdef PREDICTOR_DEBUG
  else {
    dom.getCommunicator()->fprintf(stderr, "SUCCESS: Skipped wall distance calculation (update = %d)!\n",update);
  }
  dom.getCommunicator()->fprintf(stderr, "Times: tn = %e, tnm1 = %e, tnm2 = %e\n",tPredictors[0],tPredictors[1],tPredictors[2]);
  dom.getCommunicator()->fprintf(stderr, "    tlastinit = %e, tlastinit2 = %e\n\n",tnm1,tnm2);
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
int ReinitializeDistanceToWall<dimLS>::UpdatePredictorsCheckTol(double t)
{

  int update = 0;

  if (tPredictors[1] < 0)
    update = 3;
  else if (tPredictors[2] < 0)
    update = 2;
  else {
    // compute distance predictions

  #ifdef PREDICTOR_DEBUG
    dom.getCommunicator()->fprintf(stderr,"\nUpdating predictors...\n");
  #endif

    double dtnm2, dtnm1, dtn;
    dtnm2 = tPredictors[1]-tPredictors[2];
    dtnm1 = tPredictors[0]-tPredictors[1];
    dtn = t-tPredictors[0];

    // dom.getCommunicator()->fprintf(stderr, "dtn = %e, dtnm1 = %e, dtnm2 = %e\n",dtn,dtnm1,dtnm2);

    double tolmax = 1.0e-2, tolmin = 1.0e-9, tolmin2 = 1.0e-3;
    double maxdd = -1.0;
    double delta, d2wp, deltatmp, deltam;
    int i;

  #ifdef PREDICTOR_DEBUG
    double maxdel = -1.0, meandel = 0.0, maxd2w = -1.0, mind2w = 1.0e16;
    int nvi = 0;
  #endif

  #pragma omp parallel for
    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
      for (int j = 0; j < nPredictors[iSub]; j++) {
        d2wp = (1.0+dtn/dtnm1*(dtn+2.0*dtnm1+dtnm2)/(dtnm1+dtnm2))*d2wPredictors[iSub][j][1]
        -(dtn/dtnm1)*(dtn+dtnm1+dtnm2)/dtnm2*d2wPredictors[iSub][j][2]
        +(dtn/dtnm2)*(dtn+dtnm1)/(dtnm1+dtnm2)*d2wPredictors[iSub][j][3];

        i = d2wPredictors[iSub][j][0];
        // deltatmp = 2.0*fabs(d2wp-d2wall(iSub)[i][0])/d2wall(iSub)[i][0];
        // if (d2wall(iSub)[i][0]>tolmin2)
        //   delta = deltatmp;
        // else delta = 0.0;
        delta = 2.0*fabs(d2wp-d2wall(iSub)[i][0])/d2wall(iSub)[i][0];

  #ifdef PREDICTOR_DEBUG
        // for viewing maximum predicted change and other values
        maxdel = max(delta,maxdel);
        // meandel += deltatmp;
        meandel += delta;
        mind2w = min(mind2w,d2wall(iSub)[i][0]);
        maxd2w = max(maxd2w,d2wall(iSub)[i][0]);
  #endif

  //       if (delta > tolmax) {
  //         update = 2;
  // #ifdef PREDICTOR_DEBUG
  //         fprintf(stderr,
  //           "Tolerance exceeded (2delta/d = %e > %e @ cpu %d, node %d\n    d2wp = %e, d2wall(tnm1) = %e\n    d2wnm1 = %e, d2wnm2 = %e, d2wnm3 = %e)\n",
  //           delta,tolmax,dom.getCommunicator()->cpuNum(),j,d2wp,d2wall(iSub)[i][0],d2wPredictors[iSub][j][1],d2wPredictors[iSub][j][2],d2wPredictors[iSub][j][3]);
  //         nvi++;
  // #else
  //         goto exitlbl;
  // #endif
  //       }
        // else {
          deltam = fabs(d2wp-d2wPredictors[iSub][j][1])/d2wPredictors[iSub][j][1];
          maxdd  = max(deltam,maxdd);

          d2wPredictors[iSub][j][3] = d2wPredictors[iSub][j][2];
          d2wPredictors[iSub][j][2] = d2wPredictors[iSub][j][1];
          d2wPredictors[iSub][j][1] = d2wp;
        // }
      }
    }

    // exitlbl:
    // dom.getCommunicator()->globalMax(1,&update);

  #ifdef PREDICTOR_DEBUG
    // max and mean error outputs to view (for testing)
    dom.getCommunicator()->globalMax(1,&maxdel);
    dom.getCommunicator()->fprintf(stderr,"Max delta (2deld/d) = %e vs tol = %e)\n",
      maxdel,tolmax);
    dom.getCommunicator()->globalMax(1,&maxd2w);
    dom.getCommunicator()->globalMin(1,&mind2w);
    dom.getCommunicator()->fprintf(stderr,"Max predictor d2w = %e, min = %e\n",
      maxd2w,mind2w);

    int nPredTot = 0;
    for (int iSub = 0; iSub < dom.getNumLocSub(); iSub++)
      nPredTot += nPredictors[iSub];
    dom.getCommunicator()->globalSum(1,&meandel);
    dom.getCommunicator()->globalSum(1,&nPredTot);
    meandel /= nPredTot;
    dom.getCommunicator()->fprintf(stderr,"Mean (2delta/d) = %e (N predictors = %d)\n\n",meandel,nPredTot);
  #endif

    if (meandel > tolmax) {
      update = 2;
    }

    if (update < 2) {
      dom.getCommunicator()->globalMax(1,&maxdd);
      if (maxdd < tolmin) {
  #ifdef PREDICTOR_DEBUG
        dom.getCommunicator()->fprintf(stderr,"max relative (dpn-dpnm1) = %e vs tol = %e)\n",
          maxdd,tolmin);
  #endif
        update = 1;
      }
    }
  //   else {
  // #ifdef PREDICTOR_DEBUG
  //     dom.getCommunicator()->globalSum(1,&nvi);
  //     dom.getCommunicator()->fprintf(stderr,"Total N = %d violating predictors\n",nvi);
  // #endif
  //   }

  }

  // update predictor times current timestep
  tPredictors[2] = tPredictors[1];
  tPredictors[1] = tPredictors[0];
  tPredictors[0] = t;

  dom.getCommunicator()->globalMax(1,&update);
  return update;
}

//------------------------------------------------------------------------------

// sjg, 04/2017: reinitialize predictors
template <int dimLS>
void ReinitializeDistanceToWall<dimLS>::ReinitializePredictors(int update,
                                                               double t,
                                                               DistLevelSetStructure *LSS,
                                                               DistSVec<double, 3> &X,
                                                               DistGeoState &distGeoState)
{
  tPredictors[1] = tnm1;
  tPredictors[2] = tnm2;

  int nSub = dom.getNumLocSub();

  // update distance at current predictors exactly
  if (update == 1) {

  #ifdef PREDICTOR_DEBUG
    dom.getCommunicator()->fprintf(stderr, "Reinitializing due to flat ...\n");
  #endif

    // this could cause problems if the intersections change since predictor initialization,
    // but in case of flat, shouldn't move much
    int *nPredgbg;
    nPredgbg = new int[nSub];
    sortedNodes = -1;
    tag = -1;
    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
      nSortedNodes[iSub] = 0;
      firstCheckedNode[iSub] = 0;
    }
    dom.pseudoFastMarchingMethod<1>(tag, X, d2wall, 1, 0,
      sortedNodes, nSortedNodes, firstCheckedNode, nPredgbg, LSS);

    for (int iSub = 0; iSub < nSub; ++iSub) {
      for (int j = 0; j < nPredictors[iSub]; j++) {
        int i = d2wPredictors[iSub][j][0];
        d2wPredictors[iSub][j][1] = d2wall(iSub)[i][0];
        d2wPredictors[iSub][j][2] = d2wallnm1(iSub)[i][0];
        d2wPredictors[iSub][j][3] = d2wallnm2(iSub)[i][0];
      }
    }

    delete[] nPredgbg;
  }

  // create new predictors and populate
  else if (update == 2) {

  #ifdef PREDICTOR_DEBUG
    dom.getCommunicator()->fprintf(stderr, "Reinitializing predictors...");
  #endif

    int row;
#pragma omp parallel for
    for (int iSub = 0; iSub < nSub; ++iSub) {
      row = 0;
      d2wPredictors[iSub] = new double *[nPredictors[iSub]];

      for(int i = 0; i < tag(iSub).len; i++) {
        if (tag(iSub)[i]==1) {
          d2wPredictors[iSub][row] = new double[4];
          d2wPredictors[iSub][row][0] = (double) i;
          d2wPredictors[iSub][row][1] = d2wall(iSub)[i][0];
          d2wPredictors[iSub][row][2] = d2wallnm1(iSub)[i][0];
          d2wPredictors[iSub][row][3] = d2wallnm2(iSub)[i][0];
          row++;
        }
      }

      // if (row!=nPredictors[iSub]) {
      //   fprintf(stderr, "Problem reinitializing wall distance predictors!\n");
      //   fprintf(stderr, "row = %d != nPredictors[iSub] = %d (cpu %d)\n",row,nPredictors[iSub],dom.getCommunicator()->cpuNum());
      // }
      assert (row==nPredictors[iSub]);
    }

    predictorActive = 1;

#if 0
    int printnPred = 0;
    for (int iSub = 0; iSub < nSub; iSub++)
      printnPred += nPredictors[iSub];
    dom.getCommunicator()->globalSum(1, &printnPred);
    dom.getCommunicator()->fprintf(stderr, "Initialized total of N predictors = %d\n", printnPred);
#endif
  }
  else {
      fprintf(stderr, " *** Error ***, Wall distance predictor update case is invalid\n");
      exit(1);
  }

  // set predictor times corresponding to initializations
  tPredictors[2] = tnm2;
  tPredictors[1] = tnm1;
}

//------------------------------------------------------------------------------

// The following is an adaptation of the Fast Marching Method to Embedded Turbulent computation.
// Adam 2012.09
template <int dimLS>
int ReinitializeDistanceToWall<dimLS>::PseudoFastMarchingMethod(
    DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState, int iterativeLevel, int *nPredLoc)
{
  sortedNodes = -1;
  int nSub = dom.getNumLocSub();
  if (iterativeLevel == 0)
  {
    d2wall = 1.0e10;
    tag = -1;
  }
  int isDone = 0;

  int level = iterativeLevel; // Level 0 (inActive nodes) and 1 (Embedded surface neighbors)
  while (isDone == 0)
  { // Tag and update d at every level
    dom.pseudoFastMarchingMethod<1>(tag, X, d2wall, level, iterativeLevel, sortedNodes, nSortedNodes, firstCheckedNode, nPredLoc, &LSS);
    // I don't think it is a good idea to OMP parallelize this loop. nSub should be small, though!
    isDone = 1;
    for (int iSub = 0; iSub < nSub; ++iSub)
    {
      if (nSortedNodes[iSub] != tag(iSub).len)
      {
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
void ReinitializeDistanceToWall<dimLS>::IterativeMethodUpdate(int max_level,
                                                              DistLevelSetStructure &LSS,
                                                              DistSVec<double, 3> &X,
                                                              DistGeoState &distGeoState)
{
  if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::HYBRID &&
      iod.eqs.tc.tm.d2wall.iterativelvl > 1)
    max_level = min(iod.eqs.tc.tm.d2wall.iterativelvl, max_level);

  int nSub = dom.getNumLocSub();
  bool printwarning = false;
  double maxres = -FLT_MAX;
  int maxreslvl = 1;

  int *nPredgbg;
  nPredgbg = new int[nSub];

  // new implementation !! -----------------------------------------------------

  iod.eqs.tc.tm.d2wall.eps = 1.0e-2; // correct for new norm definition

  DistSVec<double, 1> resnm1(dom.getNodeDistInfo());
  double res;
  int it;

  int nrwband = 2;
  // int nrwband = max_level-1;

  int *nSortedNodesLp, *firstCheckedNodeLp;
  nSortedNodesLp = new int[nSub];
  firstCheckedNodeLp = new int[nSub];
  DistVec<int> tagLp(dom.getNodeDistInfo());
  DistVec<int> sortedNodesLp(dom.getNodeDistInfo());

  // Propagate information outwards
  // for (int jlvl = 1; jlvl < max_level; jlvl+=nrwband) {
  for (int jlvl = 1; jlvl < max_level; jlvl++) {
    res = 1e10;
    it = 1;

    // must alter tags to allow for re-updating distance
    dom.pseudoFastMarchingMethod<1>(tag, X, d2wall, jlvl, jlvl, sortedNodes, nSortedNodes, firstCheckedNode, nPredgbg, &LSS);

    while (res > iod.eqs.tc.tm.d2wall.eps && it < iod.eqs.tc.tm.d2wall.maxIts) {
      resnm1 = d2wall;

      for (int iSub = 0; iSub < nSub; ++iSub) {
        nSortedNodesLp[iSub] = nSortedNodes[iSub];
        firstCheckedNodeLp[iSub] = firstCheckedNode[iSub];
      }
      tagLp = tag;
      sortedNodesLp = sortedNodes;

      int maxlvl = min(max_level,jlvl+nrwband);
      for (int ilvl = jlvl+1; ilvl <= maxlvl; ilvl++) {
        // Tag and update d at every level
        dom.pseudoFastMarchingMethod<1>(tagLp, X, d2wall, ilvl, jlvl, sortedNodesLp, nSortedNodesLp, firstCheckedNodeLp, nPredgbg, &LSS);
      }

      // update residual and check for convergence
      res = -1.0;
      for (int iSub = 0; iSub < nSub; ++iSub) {
        for(int i = 0; i < tag(iSub).len; i++) {
          if ((tagLp(iSub)[i]>jlvl) && (tagLp(iSub)[i]<=jlvl+nrwband))
            res = max(fabs(resnm1(iSub)[i][0]-d2wall(iSub)[i][0])/resnm1(iSub)[i][0],res);
        }
      }
      dom.getCommunicator()->globalMax(1, &res);
      it++;

      // dom.getCommunicator()->fprintf(stderr, "Residual = %e @ iteration %d\n", res, it-1);
    }

    // update tags
    tag = tagLp;

    // if (res > iod.eqs.tc.tm.d2wall.eps)
    // {
      printwarning = true;
      if (res > maxres)
      {
        maxres = res;
        maxreslvl = jlvl+1;
      }
    // }

    // dom.getCommunicator()->fprintf(stderr,
    //   "Wall distance performed %d iterations at level %d of %d (min res = %e)\n\n", it-1, jlvl+1, max_level, res);
  }
  delete[] nPredgbg;
  delete[] nSortedNodesLp;
  delete[] firstCheckedNodeLp;

  if (printwarning)
    dom.getCommunicator()->fprintf(stderr,
      "*** Warning: Iterative distance to wall computation (Max residual: %e at level: %d, target: %e)\n",
      maxres, maxreslvl, iod.eqs.tc.tm.d2wall.eps);

  // /////////// OLD IMPLEMENTATION

  // // initial residuals from first pass
  // double *res0 = new double[max_level-1];
  // std::memset(res0, 0.0, sizeof(double)*(max_level-1));
  // for (int iSub = 0; iSub < nSub; ++iSub) {
  //   for(int i = 0; i < tag(iSub).len; i++) {
  //     if(tag(iSub)[i]>1)	res0[tag(iSub)[i]-2] += d2wall(iSub)[i][0]*d2wall(iSub)[i][0];
  //   }
  // }

  // int *nSortedNodesLp, *firstCheckedNodeLp;
  // nSortedNodesLp = new int[nSub];
  // firstCheckedNodeLp = new int[nSub];
  // DistVec<int> tagLp(dom.getNodeDistInfo());
  // DistVec<int> sortedNodesLp(dom.getNodeDistInfo());

  // // Propagate information outwards
  // for (int ilvl = 1; ilvl < max_level; ++ilvl) {

  //   // double resn = res0[ilvl-1];
  //   dom.getCommunicator()->globalSum(1,res0+ilvl-1);
  //   double resn = sqrt(res0[ilvl-1]);
  //   double resnm1;
  //   double res = 1.0e10;
  //   int it = 1;

  //   // TO DO: replace calls to dom.pseudoFMM with a new function in domain which does
  //   // essentially the same by allows for update of only narrowband levels (doesn't even
  //   // need to update sortednodes list either?)

  //   // must alter tags to allow for re-updating distance
  //   dom.pseudoFastMarchingMethod<1>(tag, X, d2wall, ilvl, ilvl, sortedNodes, nSortedNodes, firstCheckedNode, nPredgbg, &LSS);

  //   while (res > iod.eqs.tc.tm.d2wall.eps && it < iod.eqs.tc.tm.d2wall.maxIts) {
  //     resnm1 = resn;

  //     for (int iSub = 0; iSub < nSub; ++iSub) {
  //       nSortedNodesLp[iSub] = nSortedNodes[iSub];
  //       firstCheckedNodeLp[iSub] = firstCheckedNode[iSub];
  //     }
  //     sortedNodesLp = sortedNodes;
  //     tagLp = tag;

  //     // calculate new distances in levels +1 and +2 (narrow band) and check conv.
  //     dom.pseudoFastMarchingMethod<1>(tagLp, X, d2wall, ilvl+1, ilvl, sortedNodesLp, nSortedNodesLp, firstCheckedNodeLp, nPredgbg, &LSS);
  //     dom.pseudoFastMarchingMethod<1>(tagLp, X, d2wall, ilvl+2, ilvl, sortedNodesLp, nSortedNodesLp, firstCheckedNodeLp, nPredgbg, &LSS);

  //     // update residual and check for convergence
  //     resn = 0.0;
  //     for (int iSub = 0; iSub < nSub; ++iSub) {
  //       for(int i = 0; i < tagLp(iSub).len; i++) {
  //         if(tagLp(iSub)[i]==ilvl+1) resn += d2wall(iSub)[i][0]*d2wall(iSub)[i][0];
  //       }
  //     }
  //     dom.getCommunicator()->globalSum(1, &resn);
  //     resn = sqrt(resn);
  //     res = fabs((resn-resnm1)/resnm1);
  //     it++;

  //     dom.getCommunicator()->fprintf(stderr, "Residual = %e @ iteration %d\n", res, it-1);
  //   }

  //   // update tags
  //   tag = tagLp;

  //   if (res > iod.eqs.tc.tm.d2wall.eps)
  //   {
  //     printwarning = true;
  //     if (res > maxres)
  //     {
  //       maxres = res;
  //       maxreslvl = ilvl+1;
  //     }
  //   }
  //   dom.getCommunicator()->fprintf(stderr,
  //     "Wall distance performed %d iterations at level %d of %d (min res = %e)\n\n", it-1, ilvl+1, max_level,res);
  // }
  // delete[] nPredgbg;
  // delete[] res0;
  // delete[] nSortedNodesLp;
  // delete[] firstCheckedNodeLp;

  // if (printwarning)
  //   dom.getCommunicator()->fprintf(stderr,
  //     "*** Warning: Iterative distance to wall computation (Max residual: %e at level: %d, target: %e)\n",
  //     maxres, maxreslvl, iod.eqs.tc.tm.d2wall.eps);
}

//------------------------------------------------------------------------------

template <int dimLS>
void ReinitializeDistanceToWall<dimLS>::InitializeWallFunction(DistLevelSetStructure &LSS,
                                                              DistSVec<double, 3> &X,
                                                              DistGeoState &distGeoState,
                                                              int *nPredLoc, double t)
{
  sortedNodes = -1;
  tag = -1;

  // Fill with initial guess
#if 1
  d2wall = 1e10;
#else
  if (t>0.0) {
  #pragma omp parallel for
    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub)
    {
      for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i)
        d2wall(iSub)[i][0] = distGeoState(iSub).getDistanceToWall()[i]>0.0?distGeoState(iSub).getDistanceToWall()[i]:1e10;
    }
  }
  else d2wall = 1e10;
#endif

  int level = 0; // Level 0 (inActive nodes) and 1 (Embedded surface neighbors)
  while (level < 2)
  { // Tag levels 0 (inactive) and 1 (embedded surface neighbors) for FIM initialization
    dom.pseudoFastMarchingMethod<1>(tag, X, d2wall, level, 0, sortedNodes, nSortedNodes, firstCheckedNode, nPredLoc, &LSS);
    ++level;
  }

  return;
}

//------------------------------------------------------------------------------

template <int dimLS>
void ReinitializeDistanceToWall<dimLS>::GetLevelsFromInterfaceAndMarchForward(DistLevelSetStructure &LSS,
                                                                              DistSVec<double, 3> &X,
                                                                              DistGeoState &distGeoState)
{
  int max_level = 1;

  // int min_level = 0;
  // int level = 1;
  int min_level = -1;
  int level = 2;  // levels <= 1 tagged in initialization already

  dummyPhi = 1.0;

  int counter = 0;
  // Tag every level
  // while (min_level <= 0)
  while (min_level < 0)
  {
    dom.TagInterfaceNodes(0, tag, dummyPhi, level, &LSS);
    min_level = 1;
    counter++;

    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub)
    {
      // for (int i = 0; i < done(iSub).len; ++i)
      for (int i = 0; i < d2wall(iSub).len; ++i)
      {
        min_level = min(min_level, tag(iSub)[i]);
        if (min_level < 0)
          goto exitlbl;
        max_level = max(max_level, tag(iSub)[i]);
      }
    }
    exitlbl:
    dom.getCommunicator()->globalMin(1, &min_level);
    // dom.getCommunicator()->fprintf(stderr, "min_level = %d (count = %d)\n",min_level,counter);
    ++level;
  }
  dom.getCommunicator()->globalMax(1, &max_level);
  // dom.getCommunicator()->fprintf(stderr, "min_level = %d, max_level = %d\n",min_level,max_level);

  if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::HYBRID &&
      iod.eqs.tc.tm.d2wall.iterativelvl > 1)
    max_level = min(iod.eqs.tc.tm.d2wall.iterativelvl, max_level);

  // Propagate information outwards
  MultiFluidData::CopyCloseNodes copy = MultiFluidData::FALSE;
  bool printwarning = false;
  double maxres = -FLT_MAX;
  int maxreslvl = 1;

  for (int ilvl = 2; ilvl <= max_level; ++ilvl) {
    double res = 1.0;
    double resn = 1.0;
    double resnm1 = 1.0;

    int it = 0;
    while (res > iod.eqs.tc.tm.d2wall.eps && it < iod.eqs.tc.tm.d2wall.maxIts)
    {
      resnm1 = resn;
      dom.computeDistanceLevelNodes(1, tag, ilvl, X, d2wall, resn, dummyPhi, copy);
      dom.getCommunicator()->globalMax(1, &resn);
      res = fabs((resn - resnm1) / (resn + resnm1));
      it++;

      // dom.getCommunicator()->fprintf(stderr, "Residual = %e @ iteration %d\n", res, it-1);
    }

    if (res > iod.eqs.tc.tm.d2wall.eps)
    {
      printwarning = true;
      if (res > maxres)
      {
        maxres = res;
        maxreslvl = ilvl;
      }
    }

    // dom.getCommunicator()->fprintf(stderr, "Wall distance performed maximum %d iterations at level %d of %d\n\n", --it, ilvl, max_level);
  }

  if (printwarning)
    dom.getCommunicator()->fprintf(stderr,
                                    "*** Warning: Iterative distance to wall computation (Max residual: %e at level: %d, target: %e)\n",
                                    maxres, maxreslvl, iod.eqs.tc.tm.d2wall.eps);
}

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
  int *nPredgbg;
  nPredgbg = new int[dom.getNumLocSub()];

  // // compare to iterative
  // int iterativeTemp = iod.eqs.tc.tm.d2wall.maxIts;
  // iod.eqs.tc.tm.d2wall.maxIts = 100;
  // dom.getCommunicator()->fprintf(stderr,
  //   "Comparing specified wall distance computation to iterative method with 100 maximum iterations.\n");
  // DistanceToClosestPointOnMovingStructure(LSS, X, distGeoState, nPredgbg);
  // GetLevelsFromInterfaceAndMarchForward(LSS, X, distGeoState);
  // iod.eqs.tc.tm.d2wall.maxIts = iterativeTemp;

  // compare to noniterative
  dom.getCommunicator()->fprintf(stderr,
    "Comparing specified wall distance computation to noniterative method.\n");
  PseudoFastMarchingMethod(LSS,X,distGeoState,0,nPredgbg);

  delete[] nPredgbg;
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
    // for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i)
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