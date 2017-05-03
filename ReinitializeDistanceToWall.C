#include <ReinitializeDistanceToWall.h>
#include <LevelSet/LevelSetStructure.h>
#include <Domain.h>
#include <DistVector.h>

// sjg, 02/2017: flags for testing and error calculations
   // #define DELTA_CHECK
   #define ERROR_CHECK
   // #define PRINT_VERB
  //  #define PRINT_INTERSECT

//------------------------------------------------------------------------------

template <int dimLS>
ReinitializeDistanceToWall<dimLS>::ReinitializeDistanceToWall(IoData &ioData, Domain &domain)
    : iod(ioData), dom(domain), done(domain.getNodeDistInfo()), d2wall(domain.getNodeDistInfo()), tag(domain.getNodeDistInfo()), dummyPhi(domain.getNodeDistInfo()), sortedNodes(domain.getNodeDistInfo()), predictorActive(0),
    d2wallnm1(domain.getNodeDistInfo()), d2wallnm2(domain.getNodeDistInfo()),tnm1(-1.0),tnm2(-1.0)
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


  dom.getCommunicator()->fprintf(stderr, "In cleanup function!\n");
  dom.getCommunicator()->barrier();

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

  // update predictors and check tolerances, receiving global maximum
  int update = UpdatePredictorsCheckTol(t);

  // predicted distance is zero, so correct predictor to exact values
  if (update == 1) {
    ReinitializePredictors(update, t, &LSS, X, distGeoState);

    d2wallnm2 = d2wallnm1;
    d2wallnm1 = d2wall;
    tnm2 = tnm1;
    tnm1 = t;
  }

  // compute distance in full domain
  else if (update > 1) {
  // if (update > 0) {

    dom.getCommunicator()->fprintf(stderr, "Performing a full domain wall distance update (update = %d)!\n",update);
    dom.getCommunicator()->barrier();



    // free memory from previous predictors if reinitializing
    // if (predictorActive && update == 2) {
    if (predictorActive) {
      for (int iSub = 0; iSub<dom.getNumLocSub(); iSub++) {
        for (int j = 0; j<nPredictors[iSub]; j++) {
          delete[] d2wPredictors[iSub][j];
        }
        delete [] d2wPredictors[iSub];
      }

      dom.getCommunicator()->fprintf(stderr, "Successfully deleted previous predictors!\n");
      dom.getCommunicator()->barrier();
    }


    // update d2wall exactly
    if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::ITERATIVE) {
      DistanceToClosestPointOnMovingStructure(LSS, X, distGeoState);
      GetLevelsFromInterfaceAndMarchForward(LSS, X, distGeoState);
    }
    else if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::NONITERATIVE) {


      dom.getCommunicator()->barrier();
      dom.getCommunicator()->fprintf(stderr, "Ready to start FMM!\n");
      dom.getCommunicator()->barrier();


      PseudoFastMarchingMethod(LSS, X, distGeoState, 0);
    }
    else if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::HYBRID) {
      int iterativeLevel = 0;

      if (iod.eqs.tc.tm.d2wall.iterativelvl > 1) {
        DistanceToClosestPointOnMovingStructure(LSS, X, distGeoState);
        GetLevelsFromInterfaceAndMarchForward(LSS, X, distGeoState);
        iterativeLevel = iod.eqs.tc.tm.d2wall.iterativelvl;
      }
      PseudoFastMarchingMethod(LSS, X, distGeoState, iterativeLevel);
    }
    else {
      fprintf(stderr, " *** Error ***, Unknown wall distance method\n");
      exit(1);
    }


    dom.getCommunicator()->fprintf(stderr, "Successfully calculated d2wall!\n");
    dom.getCommunicator()->barrier();

    // reinitialize predictors
    if (update == 2) {
    // if (update < 3) {
      ReinitializePredictors(update, t, &LSS, X, distGeoState);


      dom.getCommunicator()->fprintf(stderr, "Successfully initialized predictors!\n");
      dom.getCommunicator()->barrier();
    }

    // store previous exact updates
    d2wallnm2 = d2wallnm1;
    d2wallnm1 = d2wall;
    tnm2 = tnm1;
    tnm1 = t;

    dom.getCommunicator()->fprintf(stderr, "Successfully stored prev. d2wall!\n");
    dom.getCommunicator()->barrier();


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
    PrintIntersectedValues(&LSS, X);
#endif

  }



  dom.getCommunicator()->fprintf(stderr, "Ready to compute exact errors!\n");
  dom.getCommunicator()->barrier();


// sjg, 02/2017: wall distance error for testing
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


  dom.getCommunicator()->fprintf(stderr,"Updating times...\n");
  dom.getCommunicator()->barrier();

  int update = 0;

  // update to current timestep and compute distance predictions
  if (tPredictors[1] < 0) {
    tPredictors[2] = tPredictors[1];
    tPredictors[1] = tPredictors[0];
    tPredictors[0] = t;

    update = 3;
  }
  else if (tPredictors[2] < 0) {
    tPredictors[2] = tPredictors[1];
    tPredictors[1] = tPredictors[0];
    tPredictors[0] = t;

    update = 2;
  }
  else {


    dom.getCommunicator()->fprintf(stderr,"Updating predictors...\n");
    dom.getCommunicator()->barrier();


    double dtnm2, dtnm1, dtn;
    dtnm2 = tPredictors[1]-tPredictors[2];
    dtnm1 = tPredictors[0]-tPredictors[1];
    dtn = t-tPredictors[0];

    tPredictors[2] = tPredictors[1];
    tPredictors[1] = tPredictors[0];
    tPredictors[0] = t;

    double tolmax = 1.0e-2;
    double tolmin = 1.0e-6;
    double maxdd = -1.0;

  #pragma omp parallel for
    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
      for (int j = 0; j < nPredictors[iSub]; j++) {
        double d2wp = (1.0+dtn/dtnm1*(dtn+2.0*dtnm1+dtnm2)/(dtnm1+dtnm2))*d2wPredictors[iSub][j][1]
        -(dtn/dtnm1)*(dtn+dtnm1+dtnm2)/dtnm2*d2wPredictors[iSub][j][2]
        +(dtn/dtnm2)*(dtn+dtnm1)/(dtnm1+dtnm2)*d2wPredictors[iSub][j][3];

        d2wPredictors[iSub][j][3] = d2wPredictors[iSub][j][2];
        d2wPredictors[iSub][j][2] = d2wPredictors[iSub][j][1];
        d2wPredictors[iSub][j][1] = d2wp;

        int i = d2wPredictors[iSub][j][0];
        double delta = (2.0*(d2wp-d2wall(iSub)[i][0])/d2wall(iSub)[i][0]);

        if (delta > tolmax) {
          // dom.getCommunicator()->fprintf(stderr,
          //   "Tolerance exceeded (2delta/d = %e > %e)\n",delta,tolmax);
          fprintf(stderr,"Tolerance exceeded (2delta/d = %e > %e, d2wp = %e, d2wall = %e\n",
            delta,tolmax,d2wp,d2wall(iSub)[i][0]);
          fprintf(stderr,"d2wnm1 = %e, d2wnm2 = %e\n",
            d2wPredictors[iSub][j][2],d2wPredictors[iSub][j][3]);
          update = 2;
          goto exitlbl;
        }

        double deltam = fabs(d2wPredictors[iSub][j][1]-d2wPredictors[iSub][j][2])/d2wPredictors[iSub][j][2];
        maxdd  = max(deltam,maxdd);
      }
    }

    exitlbl:
    dom.getCommunicator()->globalMax(1,&update);

    if (update < 2) {
      dom.getCommunicator()->globalMax(1,&maxdd);
      if (maxdd < tolmin) {
        fprintf(stderr,"max relative (dpn-dpnm1) = %e < %e)\n",
          maxdd,tolmin);
        // dom.getCommunicator()->fprintf(stderr,"max relative (dpn-dpnm1) = %e < %e)\n",
        //   maxdd,tolmin);
        update = 1;
      }
    }
  }

  dom.getCommunicator()->globalMax(1,&update);
  dom.getCommunicator()->fprintf(stderr,"Returning update = %d\n",update);
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

  tPredictors[2] = tnm2;
  tPredictors[1] = tnm1;
  int nSub = dom.getNumLocSub();

  // update distance at current predictors exactly
  if (update == 1) {

  dom.getCommunicator()->fprintf(stderr, "Reinitializing due to flat ...\n");
  dom.getCommunicator()->barrier();

  int *nPredictorsgbg;
  nPredictorsgbg = new int[nSub];

#pragma omp parallel for
    for (int iSub = 0; iSub < nSub; ++iSub) {
      for (int j = 0; j < nPredictors[iSub]; j++) {
        int i = d2wPredictors[iSub][j][0];
        d2wPredictors[iSub][j][2] = d2wallnm1(iSub)[i][0];
        d2wPredictors[iSub][j][3] = d2wallnm2(iSub)[i][0];

/*** COULD SPEED UP if we instead only initialize first layer here */
/**** ALSO problem because nPredictors could change when we call FMM again */
        sortedNodes = -1;
        tag = -1;
        int level = 1;
        int iterativelevel = 1;
        dom.pseudoFastMarchingMethod<1>(tag, X, d2wall, level, iterativelevel,
          sortedNodes, nSortedNodes, firstCheckedNode, nPredictorsgbg, LSS);

        d2wPredictors[iSub][j][1] = d2wall(iSub)[i][0];
      }
    }

    delete[] nPredictorsgbg;
  }

  // create new predictors and populate
  else if (update == 2) {


  dom.getCommunicator()->fprintf(stderr, "Reinitializing predictors...\n");
  dom.getCommunicator()->barrier();


#pragma omp parallel for
    int row;
    double **vtag;
    vtag = new double*[nSub];
    for (int iSub = 0; iSub < nSub; ++iSub) {
      row = 0;
      vtag[iSub] = new double[distGeoState(iSub).getDistanceToWall().size()];
      std::memset(vtag[iSub], 0, sizeof(double)*distGeoState(iSub).getDistanceToWall().size());

      d2wPredictors[iSub] = new double *[nPredictors[iSub]];
      LevelSetStructure *locLSS = &((*LSS)(iSub));

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
            d2wPredictors[iSub][row] = new double[4];
            d2wPredictors[iSub][row][0] = (double) i;
            d2wPredictors[iSub][row][1] = d2wall(iSub)[i][0];
            d2wPredictors[iSub][row][2] = d2wallnm1(iSub)[i][0];
            d2wPredictors[iSub][row][3] = d2wallnm2(iSub)[i][0];
            row++;
            vtag[iSub][i] = 1;
          }
          if (jActive && !vtag[iSub][j]) {
            d2wPredictors[iSub][row] = new double[4];
            d2wPredictors[iSub][row][0] = (double) j;
            d2wPredictors[iSub][row][1] = d2wall(iSub)[j][0];
            d2wPredictors[iSub][row][2] = d2wallnm1(iSub)[j][0];
            d2wPredictors[iSub][row][3] = d2wallnm2(iSub)[j][0];
            row++;
            vtag[iSub][j] = 1;
          }
        }
      }
      // DEBUG CHECK
      if (row != nPredictors[iSub])
        fprintf(stderr, "Problem: nPredictors is incorrect or not all first layer nodes are captured!\n");

      delete[] vtag[iSub];
    }

    predictorActive = 1;
    delete[] vtag;

  }
  else {
      fprintf(stderr, " *** Error ***, Wall distance predictor update case is invalid\n");
      exit(1);
  }
}

//------------------------------------------------------------------------------

template <int dimLS>
void ReinitializeDistanceToWall<dimLS>::DistanceToClosestPointOnMovingStructure(DistLevelSetStructure &LSS,
                                                                                DistSVec<double, 3> &X,
                                                                                DistGeoState &distGeoState)
{
  done = false;
  tag = 0;

#pragma omp parallel for
  for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub)
  {

// Fill with initial guess
#if 1
    d2wall = 1e10;
#else
    for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i)
      d2wall(iSub)[i][0] = distGeoState(iSub).getDistanceToWall()[i];
#endif

    InitializeWallFunction(*dom.getSubDomain()[iSub], LSS(iSub), done(iSub), X(iSub), d2wall(iSub), tag(iSub));
    dom.getSubDomain()[iSub]->sndData(*dom.getVolPat(), d2wall(iSub).data());
  }

  dom.getVolPat()->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub)
    dom.getSubDomain()[iSub]->minRcvData(*dom.getVolPat(), d2wall(iSub).data());

  return;

// TAKEN FROM PSEUDOFMM, REUSE CODE TO INITIALIZE!
  // sortedNodes = -1;
  // d2wall = 1.0e10;
  // tag = -1;

  // int level = iterativeLevel; // Level 0 (inActive nodes) and 1 (Embedded surface neighbors)
  // while (level < 2)
  // { // Tag every level
  //   dom.pseudoFastMarchingMethod<1>(tag, X, d2wall, level, iterativeLevel, sortedNodes, nSortedNodes, firstCheckedNode, &LSS);
  //   ++level;
  // }
  // dom.getCommunicator()->globalMax(1, &level);

}

//------------------------------------------------------------------------------

template <int dimLS>
void ReinitializeDistanceToWall<dimLS>::InitializeWallFunction(SubDomain &subD,
                                                               LevelSetStructure &LSS,
                                                               Vec<bool> &done, SVec<double, 3> &X,
                                                               SVec<double, 1> &d2w, Vec<int> &tag)
{
  int(*ptrEdge)[2] = subD.getEdges().getPtr();

  for (int l = 0; l < subD.getEdges().size(); ++l)
  {
    if (LSS.edgeIntersectsStructure(0, l))
    {
      int i = ptrEdge[l][0];
      int j = ptrEdge[l][1];

      done[i] = true;
      tag[i] = 1;

      LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0, l, true);
      d2w[i][0] = LSS.isPointOnSurface(X[i], resij.trNodes[0], resij.trNodes[1], resij.trNodes[2]);

      done[j] = true;
      tag[j] = 1;

      LevelSetResult resji = LSS.getLevelSetDataAtEdgeCenter(0.0, l, false);
      d2w[j][0] = LSS.isPointOnSurface(X[j], resji.trNodes[0], resji.trNodes[1], resji.trNodes[2]);
    }
  }
}

//------------------------------------------------------------------------------

template <int dimLS>
void ReinitializeDistanceToWall<dimLS>::GetLevelsFromInterfaceAndMarchForward(DistLevelSetStructure &LSS,
                                                                              DistSVec<double, 3> &X,
                                                                              DistGeoState &distGeoState)
{
  int max_level = 1;
  int min_level = 0;
  int level = 1;

  dummyPhi = 1.0;

  // Tag every level
  while (min_level <= 0)
  {
    dom.TagInterfaceNodes(0, tag, dummyPhi, level, &LSS);
    min_level = 1;

    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub)
    {
      for (int i = 0; i < done(iSub).len; ++i)
      {
        min_level = min(min_level, tag(iSub)[i]);
        max_level = max(max_level, tag(iSub)[i]);
      }
    }
    dom.getCommunicator()->globalMin(1, &min_level);
    ++level;
  }
  dom.getCommunicator()->globalMax(1, &max_level);

  if (iod.eqs.tc.tm.d2wall.type == WallDistanceMethodData::HYBRID &&
      iod.eqs.tc.tm.d2wall.iterativelvl > 1)
    max_level = min(iod.eqs.tc.tm.d2wall.iterativelvl, max_level);

  // Propagate information outwards
  MultiFluidData::CopyCloseNodes copy = MultiFluidData::FALSE;
  bool printwarning = false;
  double maxres = -FLT_MAX;
  int maxreslvl = 1;

  for (int ilvl = 2; ilvl <= max_level; ++ilvl)
  {
    double res = 1.0;
    double resn = 1.0;
    double resnm1 = 1.0;

    int it = 0;

    while (res > iod.eqs.tc.tm.d2wall.eps && it < iod.eqs.tc.tm.d2wall.maxIts)
    {
      resnm1 = resn;
      dom.computeDistanceLevelNodes(1, tag, ilvl, X, d2wall, resn, dummyPhi, copy);
      dom.getCommunicator()->globalMax(1, &resn);
      it++;
      res = fabs((resn - resnm1) / (resn + resnm1));
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

// sjg, 02/2017: debugging output number of iterations
#ifdef PRINT_VERB
    dom.getCommunicator()->fprintf(stderr, "Wall distance performed %d iterations at level %d of %d\n", --it, ilvl, max_level);
#endif
  }

  if (printwarning)
    dom.getCommunicator()->fprintf(stderr,
                                   "*** Warning: Distance to wall computation (Max residual: %e at level: %d, target: %e)\n",
                                   maxres, maxreslvl, iod.eqs.tc.tm.d2wall.eps);
}

//------------------------------------------------------------------------------

// The following is an adaptation of the Fast Marching Method to Embedded Turbulent computation.
// Adam 2012.09
template <int dimLS>
void ReinitializeDistanceToWall<dimLS>::PseudoFastMarchingMethod(
    DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState, int iterativeLevel)
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
  { // Tag every level
    dom.pseudoFastMarchingMethod<1>(tag, X, d2wall, level, iterativeLevel, sortedNodes, nSortedNodes, firstCheckedNode, nPredictors, &LSS);
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
  dom.getCommunicator()->globalMax(1, &level);

// sjg, 02/2017: wall distance print number of levels
#ifdef PRINT_VERB
  dom.getCommunicator()->fprintf(stderr, "There are %d levels\n", --level);
#endif

  return;
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

  // if(iod.eqs.tc.tm.d2wall.type ==  WallDistanceMethodData::NONITERATIVE) {
  int iterativeTemp = iod.eqs.tc.tm.d2wall.maxIts;
  iod.eqs.tc.tm.d2wall.maxIts = 100;
  dom.getCommunicator()->fprintf(stderr,
    "Comparing specified wall distance computation to iterative method with 100 maximum iterations.\n");
  DistanceToClosestPointOnMovingStructure(LSS, X, distGeoState);
  GetLevelsFromInterfaceAndMarchForward(LSS, X, distGeoState);
  dom.getCommunicator()->fprintf(stderr,
    "Computing wall distance errors.\n");
  iod.eqs.tc.tm.d2wall.maxIts = iterativeTemp;
  // }
  // else if(iod.eqs.tc.tm.d2wall.type ==  WallDistanceMethodData::ITERATIVE) {
  //   PseudoFastMarchingMethod(LSS,X,distGeoState,0);
  //   dom.getCommunicator()->fprintf(stderr,"Comparing specified wall distance computation to non-iterative method.\n");
  // }
  DistSVec<double, 1> d2wall_ref = d2wall;

#pragma omp parallel for
  for (int iSub = 0; iSub < nSub; ++iSub)
  {
    errors[iSub] = new double[4];
    errors[iSub][0] = 0.0;          // mean
    errors[iSub][1] = 0.0;          // RMS
    errors[iSub][2] = 0.0;          // max
    errors[iSub][3] = 0.0;          // distance of max
    // errors[iSub][4] = 0.0; errors[iSub][5] = 0.0;
    errorsEx[iSub] = new double[4];
    errorsEx[iSub][0] = 0.0;
    errorsEx[iSub][1] = 0.0;
    errorsEx[iSub][2] = 0.0;
    errorsEx[iSub][3] = 0.0;
    // errorsEx[iSub][4] = 0.0; errorsEx[iSub][5] = 0.0;

    nDofs[iSub] = 0;
    for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i)
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
          // errors[iSub][4] = iSub;
          // errors[iSub][5] = i;
        }
        localErrorEx = fabs(d2wall_ref(iSub)[i][0] - d2wall_comp(iSub)[i][0]) / d2wall_ref(iSub)[i][0];
        errorsEx[iSub][0] += localErrorEx;
        errorsEx[iSub][1] += localErrorEx * localErrorEx;
        if (localErrorEx > errorsEx[iSub][2])
        {
          errorsEx[iSub][2] = localErrorEx;
          errorsEx[iSub][3] = d2wall_ref(iSub)[i][0];
          // errorsEx[iSub][4] = iSub;
          // errorsEx[iSub][5] = i;
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
      // errors[0][4] = errors[iSub][4];
      // errors[0][5] = errors[iSub][5];
    }
    errorsEx[0][0] += errorsEx[iSub][0];
    errorsEx[0][1] += errorsEx[iSub][1];
    if (errorsEx[iSub][2] > errorsEx[0][2])
    {
      errorsEx[0][2] = errorsEx[iSub][2];
      errorsEx[0][3] = errorsEx[iSub][3];
      // errorsEx[0][4] = errorsEx[iSub][4];
      // errorsEx[0][5] = errorsEx[iSub][5];
    }
  }

  // Communicate across all processes to find global sums/max
  dom.getCommunicator()->globalSum(1, nDofs);
  dom.getCommunicator()->globalSum(2, errors[0]);
  dom.getCommunicator()->globalSum(2, errorsEx[0]);

  MPI_Comm comm = dom.getCommunicator()->getMPIComm();
  MPI_Allreduce(&errors[0][2], &errors[0][2], 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, comm);
  MPI_Allreduce(&errorsEx[0][2], &errorsEx[0][2], 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, comm);
  // dom.getCommunicator()->globalMax(1,errors[0]+2);
  // dom.getCommunicator()->globalMax(1,errorsEx[0]+2);

  errors[0][0] /= nDofs[0];
  errors[0][1] /= nDofs[0];
  errors[0][1] = sqrt(errors[0][1]);
  errorsEx[0][0] /= nDofs[0];
  errorsEx[0][1] /= nDofs[0];
  errorsEx[0][1] = sqrt(errorsEx[0][1]);

  dom.getCommunicator()->fprintf(stderr, "Absolute d2wall Error: %12.8e, %12.8e, %12.8e at %12.8e\n", errors[0][0], errors[0][1], errors[0][2], errors[0][3]);
  dom.getCommunicator()->fprintf(stderr, "Relative d2wall Error: %12.8e, %12.8e, %12.8e at %12.8e\n\n", errorsEx[0][0], errorsEx[0][1], errorsEx[0][2], errorsEx[0][3]);
  // dom.getCommunicator()->fprintf(stderr,"Absolute d2wall Error: %12.8e, %12.8e, %12.8e at %12.8e (iSub = %8.8d, i = %8.8d)\n",errors[0][0],errors[0][1],errors[0][2],errors[0][3],(int) errors[0][4],(int) errors[0][5]);
  // dom.getCommunicator()->fprintf(stderr,"Relative d2wall Error: %12.8e, %12.8e, %12.8e at %12.8e (iSub = %8.8d, i = %8.8d)\n",errorsEx[0][0],errorsEx[0][1],errorsEx[0][2],errorsEx[0][3],(int) errorsEx[0][4],(int) errorsEx[0][5]);

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
    // d2wChange[iSub][4] = 0.0;
    // d2wChange[iSub][5] = 0.0;

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
          // d2wChange[iSub][4] = iSub;
          // d2wChange[iSub][5] = i;
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
      // d2wChange[0][4] = d2wChange[iSub][4];
      // d2wChange[0][5] = d2wChange[iSub][5];
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
  // dom.getCommunicator()->fprintf(stderr,"Maximum absolute d2wall change since last call: %12.8e at %12.8e (iSub = %8.8d, i = %8.8d)\n",d2wChange[0][2],d2wChange[0][3],(int)d2wChange[0][4],(int)d2wChange[0][5]);

  for (int iSub = 0; iSub < nSub; iSub++)
    delete[] d2wChange[iSub];
  delete[] d2wChange;
  return;
}

//------------------------------------------------------------------------------

// sjg, 04/2017: print wall distance at intersected edges
template <int dimLS>
void ReinitializeDistanceToWall<dimLS>::PrintIntersectedValues(DistLevelSetStructure *LSS, DistSVec<double, 3> &X)
{
  dom.getCommunicator()->barrier();
  dom.getCommunicator()->fprintf(stderr, "Printing distance at intersected edges\n\n");
  dom.getCommunicator()->barrier();
  for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
    LevelSetStructure *locLSS = &((*LSS)(iSub));

    int(*ptrEdge)[2] = (*dom.getSubDomain()[iSub]).getEdges().getPtr();
    for (int l = 0; l < (*dom.getSubDomain()[iSub]).getEdges().size(); ++l)
    {
      if (locLSS->edgeIntersectsStructure(0, l))
      {
        int i = ptrEdge[l][0];
        int j = ptrEdge[l][1];
        bool iActive = locLSS->isActive(0.0,i);
        bool jActive = locLSS->isActive(0.0,j);

        if (iActive) {
          Vec3D pt = X(iSub)[i];
          fprintf(stderr, "d2w(cpu %3d, lsubdom %3d, node %5d at %12.8e,%12.8e,%12.8e) = %12.8e\n",
            dom.getCommunicator()->cpuNum(),iSub,i,pt[0],pt[1],pt[2],d2wall(iSub)[i][0]);
        }
        if (jActive) {
          Vec3D pt = X(iSub)[j];
          fprintf(stderr, "d2w(cpu %3d, lsubdom %3d, node %5d at %12.8e,%12.8e,%12.8e) = %12.8e\n",
            dom.getCommunicator()->cpuNum(),iSub,j,pt[0],pt[1],pt[2],d2wall(iSub)[j][0]);
        }
      }
    }
  }
  dom.getCommunicator()->barrier();
  dom.getCommunicator()->fprintf(stderr, "\n");
  dom.getCommunicator()->barrier();
}

//------------------------------------------------------------------------------
template class ReinitializeDistanceToWall<1>;
template class ReinitializeDistanceToWall<2>;