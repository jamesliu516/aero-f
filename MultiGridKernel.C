/* MultiGridKernel.C

 */

#include <MultiGridKernel.h>
//#define MG_TEST_LEVEL

#ifndef MG_TEST_LEVEL
static int numSmooths_pre[] = {1,2,3,3,3,3,3,3,3,3,3,3};
#else
static int numSmooths_pre[] = {0,500,0,1000,0,0,0,0,0,0};
#endif

static int numSmooths_post[] = {0,0,0,0,0,0,0,0,0,0,0};

template<class Scalar,int dim>
class MultiGridPrecMatVecProd {

  DistMat<Scalar,dim>* macroA;

  MultiGridLevel<Scalar>* multiGridLevel;

 public:

  MultiGridPrecMatVecProd(DistMat<Scalar,dim>* macroA, MultiGridLevel<Scalar>* multiGridLevel) :
                          macroA(macroA), multiGridLevel(multiGridLevel) { }

  ~MultiGridPrecMatVecProd() { }

  void apply(DistSVec<Scalar,dim>& x, DistSVec<Scalar,dim>& b) {

    b = 0.0;

    multiGridLevel->computeMatVecProd(*dynamic_cast< DistMvpMatrix<Scalar,dim>*>(macroA),
                                      x, b);
    
  }
};

template<class Scalar,int dim,class Scalar2>
class MultiGridPrecJacobiPrec {

  MultiGridSmoothingMatrix<Scalar,dim>** smoothingMatrices;
  MultiGridLevel<Scalar>* level;

 public:

  MultiGridPrecJacobiPrec(MultiGridSmoothingMatrix<Scalar,dim>** smoothingMatrices,
                          MultiGridLevel<Scalar>* level) :
                         smoothingMatrices(smoothingMatrices),
                         level(level) { }

  ~MultiGridPrecJacobiPrec() { }

  void apply(DistSVec<Scalar2,dim>& x, DistSVec<Scalar2,dim>& b) {

    b = 0.0;
    int numLocSub = x.info().numLocSub;
#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub) {
      smoothingMatrices[iSub]->smooth(x(iSub),b(iSub));
    }
    level->assemble(b);
  }
};

template<class Scalar,int dim,class Scalar2>
class MultiGridPrecRASPrec {

  MultiGridSmoothingMatrix<Scalar,dim>** smoothingMatrices;
  MultiGridLevel<Scalar>* level;

 public:

  MultiGridPrecRASPrec(MultiGridSmoothingMatrix<Scalar,dim>** smoothingMatrices,
                          MultiGridLevel<Scalar>* level) :
                         smoothingMatrices(smoothingMatrices),
                         level(level) { }

  ~MultiGridPrecRASPrec() { }

  void apply(DistSVec<Scalar2,dim>& x, DistSVec<Scalar2,dim>& b) {

    b = 0.0;
    int numLocSub = x.info().numLocSub;
#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub) {
      smoothingMatrices[iSub]->smooth(x(iSub),b(iSub));
    }
    level->assemble(b);
  }
};


template<class Scalar, int dim, class Scalar2>
MultiGridKernel<Scalar,dim,Scalar2>::MultiGridKernel(Domain *dom, DistGeoState& distGeoState, KspData& coarseSolverData, IoData& ioData,VarFcn* varFcn,bool createFineA,int num_levels, DistTimeState<dim>* ts,BcFcn* bcFcn)
    :  domain(dom), num_levels(num_levels), agglom_size(8), numLocSub(dom->getNumLocSub()), multiGridLevels(new MultiGridLevel<Scalar2>*[num_levels+1]),
    geoState(distGeoState), initialized(false), coarseSolverData(coarseSolverData),ioData(ioData), myVarFcn(varFcn),myTimeState(ts), defaultSmoother(this), bcFcn(bcFcn)
{

  ownsFineA = createFineA;

  isGeometric = false;
  myOperators = new MultiGridOperator<Scalar2,dim>*[num_levels+1];
  for (int lvl = 0; lvl < num_levels+1; ++lvl)
    myOperators[lvl] = NULL;

  mySmoother = &defaultSmoother;

  prolong_relax_factor = restrict_relax_factor = 1.0;

  beta = ioData.mg.directional_coarsening_factor;
}

template<class Scalar, int dim, class Scalar2>
void MultiGridKernel<Scalar,dim,Scalar2>::setParameters(int v1, int v2, int
finesweeps, double relax, int do_out) {

  nSmooth1 = v1;
  nSmooth2 = v2;
  relaxationFactor = relax;
  fine_sweeps = finesweeps;
  output = do_out;
}

template<class Scalar, int dim, class Scalar2>
void MultiGridKernel<Scalar,dim,Scalar2>::setOperators(SpaceOperator<dim>* spo) {

  myFluxFcn = spo->getFluxFcn();
  myVarFcn = spo->getVarFcn(); 
}

template<class Scalar, int dim, class Scalar2>
void MultiGridKernel<Scalar,dim,Scalar2>::initialize() {

  initialized = true;

  MultiGridMethod mgm = MultiGridGeometric;

  smoothingMatrices = new MultiGridSmoothingMatrix<Scalar2,dim>**[num_levels+1];

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    domain->getSubDomain()[iSub]->getEdges().updateLength(geoState.getXn()(iSub));
    domain->getSubDomain()[iSub]->makeMasterFlag(domain->getNodeDistInfo());
  }

  multiGridLevels[0] = new MultiGridLevel<Scalar2>(mgm,NULL,*domain, domain->getNodeDistInfo(), domain->getEdgeDistInfo());
  multiGridLevels[0]->copyRefinedState(domain->getNodeDistInfo(), domain->getEdgeDistInfo(),domain->getInletNodeDistInfo(), domain->getFaceDistInfo(),geoState, *domain);

  char top_file_name[256];
  for(int level = 0; level < num_levels; ++level) {
    multiGridLevels[level+1] = new MultiGridLevel<Scalar2>(mgm,multiGridLevels[level],*domain, multiGridLevels[level]->getNodeDistInfo(), multiGridLevels[level]->getEdgeDistInfo());
    multiGridLevels[level+1]->agglomerate(multiGridLevels[level]->getNodeDistInfo(),
                                          multiGridLevels[level]->getEdgeDistInfo(),
                                          multiGridLevels[level]->getIdPat(),
                                          multiGridLevels[level]->getSharedNodes(),
                                          multiGridLevels[level]->getConnectivity(),
                                          multiGridLevels[level]->getEdges(),
                                          multiGridLevels[level]->getSharedEdges(),
                                          multiGridLevels[level]->getNumSharedEdges(),
                                          *domain,dim,
                                          multiGridLevels[level]->getEdgeNormals(),
                                          multiGridLevels[level]->getCtrlVol(),
                                          NULL,beta);
    sprintf(top_file_name,"level%d",level+1);
    multiGridLevels[level+1]->writePVTUFile(top_file_name);
    sprintf(top_file_name,"agglevel%d",level+1);
    multiGridLevels[level+1]->writePVTUAgglomerationFile(top_file_name);
    domain->getCommunicator()->fprintf(stdout,"Agglomerated level %d\n", level+1);
    fflush(stdout);
  }

  macroA = new DistMat<Scalar2,dim>*[num_levels+1];

  macroValues = new DistSVec<Scalar2,dim>*[num_levels+1];
  macroValuesOld = new DistSVec<Scalar2,dim>*[num_levels+1];
  macroIrey = new DistVec<Scalar2>*[num_levels+1];
  macroValuesTmp = new DistSVec<Scalar2,dim>*[num_levels+1];
  macroR = new DistSVec<Scalar2,dim>*[num_levels+1];
  macroDX = new DistSVec<Scalar2,dim>*[num_levels+1];
  macroValuesTmpV = new DistSVec<Scalar2,dim>*[num_levels+1];
  
  for(int level = 0; level <= num_levels; ++level) {
    macroA[level] = 0;
    if (level > 0 || ownsFineA)
      macroA[level] = new DistMvpMatrix<Scalar2,dim>(domain,multiGridLevels[level]->getNodeDistInfo().subLen,
                                                  multiGridLevels[level]->getEdges(),
                                                  NULL); 

    macroValues[level] = new DistSVec<Scalar2,dim>(multiGridLevels[level]->getNodeDistInfo());
    macroValuesTmpV[level] = new DistSVec<Scalar2,dim>(multiGridLevels[level]->getNodeDistInfo());
    macroValuesOld[level] = new DistSVec<Scalar2,dim>(multiGridLevels[level]->getNodeDistInfo());
    macroValuesTmp[level] = new DistSVec<Scalar2,dim>(multiGridLevels[level]->getNodeDistInfo());
    macroR[level] = new DistSVec<Scalar2,dim>(multiGridLevels[level]->getNodeDistInfo());
    macroIrey[level] = new DistVec<Scalar2>(multiGridLevels[level]->getNodeDistInfo());
    macroDX[level] = new DistSVec<Scalar2,dim>(multiGridLevels[level]->getNodeDistInfo());
  }

  typename MultiGridSmoothingMatrix<Scalar2,dim>::SmoothingMode smoothing_mode;
  /*if (pcData.mg_smoother == PcData::MGJACOBI)
    smoothing_mode = MultiGridSmoothingMatrix<Scalar2,dim>::BlockJacobi;
  else if (pcData.mg_smoother == PcData::MGLINEJACOBI)
    smoothing_mode = MultiGridSmoothingMatrix<Scalar2,dim>::LineJacobi;
  else if (pcData.mg_smoother == PcData::MGRAS) */
    smoothing_mode = MultiGridSmoothingMatrix<Scalar2,dim>::RAS;

  for(int level = 0; level <= num_levels; ++level) {
    smoothingMatrices[level] = new MultiGridSmoothingMatrix<Scalar2,dim>*[numLocSub];
#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub)
      smoothingMatrices[level][iSub] = new
        MultiGridSmoothingMatrix<Scalar2,dim>((level == num_levels ? MultiGridSmoothingMatrix<Scalar2,dim>::RAS : smoothing_mode),
                                              iSub,
                                              multiGridLevels[level]->getNodeDistInfo().subSize(iSub),
                                              multiGridLevels[level]->getEdgeDistInfo().subSize(iSub),
                                              0,
                                              (level < num_levels ? multiGridLevels[level+1]: NULL),
                                              multiGridLevels[level]);
  }

  coarseMvps = new MultiGridPrecMatVecProd<Scalar2,dim>*[num_levels];
  coarsePrecs = new MultiGridPrecRASPrec<Scalar2,dim>*[num_levels];
  coarseSolvers = new  KspSolver<DistSVec<Scalar2,dim>, MultiGridPrecMatVecProd<Scalar2,dim>,
                                   MultiGridPrecRASPrec<Scalar2,dim>, Communicator>*[num_levels];

  for (int lvl = 1; lvl <= num_levels; ++lvl) {
    coarseMvps[lvl] = new MultiGridPrecMatVecProd<Scalar2,dim>(macroA[lvl], 
                                                        multiGridLevels[lvl]);

    coarsePrecs[lvl] = new MultiGridPrecRASPrec<Scalar2,dim>(smoothingMatrices[lvl],
                                                          multiGridLevels[lvl]);

    coarseSolvers[lvl] = new GmresSolver<DistSVec<Scalar2,dim>, MultiGridPrecMatVecProd<Scalar2,dim>,
                                   MultiGridPrecRASPrec<Scalar2,dim>, Communicator>(
                       macroValues[lvl]->info(), coarseSolverData, coarseMvps[lvl],
                       coarsePrecs[lvl], domain->getCommunicator());

  }
/*  coarseSolver->setEps(1.0e-12);
  coarseSolver->setMaxIts(100);
  coarseSolver->disableOutput();
*/
/*#pragma omp parallel for
      for(int iSub = 0; iSub < numLocSub; ++iSub) {

        for (int k = 0; k < (*macroR[0])(iSub).size(); ++k)
          tmp(iSub)[k][0] = static_cast<double>((multiGridLevels[1]->getNodeMapping())(iSub)[k]);
     }   
      char fn[32];
     
   sprintf(fn,"nodeMapping");
        domain->writeVectorToFile(fn,0,(double)0.0,tmp,
                                  &scale);
  
  
  sprintf(top_file_name,"level%i.top",0);
  multiGridLevels[0]->WriteTopFile(top_file_name);
  for(int level = 0; level < num_levels; ++level) {
    //multiGridLevels[level+1]->computeRestrictedQuantities(multiGridLevels[level]->getDistGeoState());
    sprintf(top_file_name,"level%i.top",level+1);
    multiGridLevels[level+1]->WriteTopFile(top_file_name);
  }
*/
}

template<class Scalar, int dim, class Scalar2>
MultiGridKernel<Scalar,dim,Scalar2>::~MultiGridKernel()
{
#pragma omp parallel for
  for (int level = 0; level <= num_levels; ++level) {
    delete macroValues[level];
    delete macroValuesTmp[level];
    delete macroR[level];
    delete macroDX[level];
    delete multiGridLevels[level];
  }
  delete []macroA;
  delete []macroValues;
  delete []macroValuesTmp;
  delete []macroR;
  delete []macroDX;
  delete []multiGridLevels;
}

template<class Scalar, int dim, class Scalar2>
void MultiGridKernel<Scalar,dim,Scalar2>::getData(DistMat<Scalar2,dim>& mat) {

  macroA[0] = &mat;
}

template<class Scalar, int dim, class Scalar2>
void MultiGridKernel<Scalar,dim,Scalar2>::setupAlgebraic()
{
//  char top_file_name[256];
//  sprintf(top_file_name,"level%i.top",0);
//  multiGridLevels[0]->WriteTopFile(top_file_name);
  for(int level = 0; level < num_levels; ++level) {
  //  multiGridLevels[level+1]->computeRestrictedQuantities(multiGridLevels[level]->getDistGeoState());
//    sprintf(top_file_name,"level%i.top",level+1);
//    multiGridLevels[level+1]->WriteTopFile(top_file_name);
  }

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    smoothingMatrices[0][iSub]->getData((*macroA[0])(iSub));
  }

  for(int level = 0; level < num_levels; ++level) {
    multiGridLevels[level+1]->RestrictOperator(*multiGridLevels[level],*macroA[level],
                                               *macroA[level+1]);
    if (level < num_levels) {
#pragma omp parallel for
      for(int iSub = 0; iSub < numLocSub; ++iSub) {
        smoothingMatrices[level+1][iSub]->getData((*macroA[level+1])(iSub));
      }
    }
  }
}

template<class Scalar, int dim, class Scalar2>
void MultiGridKernel<Scalar,dim,Scalar2>::
setupPreconditioner(int level) {

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    smoothingMatrices[level][iSub]->getData((*macroA[level])(iSub));
  }
}

template<class Scalar, int dim, class Scalar2>
void MultiGridKernel<Scalar,dim,Scalar2>::
setGeometric() {

  isGeometric = true;
  // In the case we are doing geometric multigrid, we need to construct operators
  for (int lvl = 1; lvl < num_levels+1; ++lvl)
    myOperators[lvl] = new MultiGridOperator<Scalar2,dim>(multiGridLevels[lvl],
                                                          ioData, myVarFcn, 
                                                          domain, bcFcn);  
}

template<class Scalar, int dim, class Scalar2>
void MultiGridKernel<Scalar,dim,Scalar2>::
setupGeometric(DistSVec<Scalar2,dim>& U) {

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    smoothingMatrices[0][iSub]->getData((*macroA[0])(iSub));
  }

  *macroValues[0] = U;
  *macroIrey[0] = *myTimeState->getInvReynolds();
//  multiGridLevels[0]->writeXpostFile("dt0",myTimeState->getDt());
  for (int lvl = 1; lvl < num_levels+1; ++lvl) {
 
    DistTimeState<dim>* fineTimeState = (lvl==1?myTimeState:&myOperators[lvl-1]->getTimeState());
    myOperators[lvl]->getTimeState().copyTimeData(fineTimeState); 
    multiGridLevels[lvl]->Restrict(*multiGridLevels[lvl-1],
                                   *macroValues[lvl-1],
                                   *macroValues[lvl]);
    multiGridLevels[lvl]->Restrict(*multiGridLevels[lvl-1],
                                   *macroIrey[lvl-1],
                                   *macroIrey[lvl]);
    multiGridLevels[lvl]->Restrict(*multiGridLevels[lvl-1],
                                   fineTimeState->getDt(),
                                   myOperators[lvl]->getTimeState().getDt());
    myVarFcn->conservativeToPrimitive(*macroValues[lvl], *macroValuesTmp[lvl]);
    myOperators[lvl]->computeJacobian(*macroValues[lvl],
                                      *macroValuesTmp[lvl],
//                                      *macroIrey[lvl],
                                      myFluxFcn,
                                      dynamic_cast< DistMvpMatrix<Scalar2,dim>&>(*macroA[lvl]));
//    multiGridLevels[lvl]->writeXpostFile("dt",myOperators[lvl]->getTimeState().getDt());
#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub) {
      smoothingMatrices[lvl][iSub]->getData((*macroA[lvl])(iSub));
    }
  }
}

template<class Scalar, int dim, class Scalar2>
void MultiGridKernel<Scalar,dim,Scalar2>::setupBcs(DistBcData<dim>* bcd) {

  for(int level = 1; level <= num_levels; ++level) {
 
    if (level == 1) 
      multiGridLevels[level]->setupBcs(*bcd, myOperators[level]->getBcData(),myOperators[level]->getBoundaryState());
    else
      multiGridLevels[level]->setupBcs(myOperators[level-1]->getBcData(), myOperators[level]->getBcData(),myOperators[level]->getBoundaryState());
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void MultiGridKernel<Scalar,dim,Scalar2>::cycleV(DistSVec<Scalar2,dim>& f, 
                                                DistSVec<Scalar2,dim>& x) {

  static int u_tag = 0;
  double scale = 1.0;
  static int cycle_counts[] = {0,0,0,0,0,0,0};
  *macroValues[0] = x;
  mySmoother->smooth(0, *macroValues[0], f,numSmooths_pre[0]);
  *macroValuesTmp[0] = f;
  *macroValuesOld[0] = 0.0;//*macroValues[0];

//  multiGridLevels[0]->writeXpostFile("res00",*macroR[0],0); 
 
  for(int level = 0; level < num_levels; ++level) {
    if (!isGeometric)
      *macroValues[level+1] = 0.0; 
    else {
      multiGridLevels[level+1]->Restrict(*multiGridLevels[level], *macroValues[level],
                                        *macroValues[level+1]);
    }

    // Restrict the residual. (I_h*d_h)
    multiGridLevels[level+1]->Restrict(*multiGridLevels[level], *macroR[level], *macroR[level+1]);
//    multiGridLevels[1]->writeXpostFile("res10",*macroR[1],4); 
    *macroValuesOld[level+1] = *macroValues[level+1];
    if (isGeometric) {
      mySmoother->applyOperator(level+1, *macroValues[level+1], *macroValuesTmp[level+1]);
    } else {

      *macroValuesTmp[level+1] = 0.0;
    }

#ifdef MG_TEST_LEVEL 
    *macroValuesTmp[level+1] = 0.0;
#else   
    //multiGridLevels[level+1]->ProjectResidual(*macroR[level+1]);
    *macroValuesTmp[level+1] += *macroR[level+1]*restrict_relax_factor;
    myOperators[level+1]->applyBCsToResidual(*macroValues[level+1], *macroValuesTmp[level+1]);
    int r;
    MPI_Comm_rank(MPI_COMM_WORLD,&r);
    if (r == 237 && level == 0) {

      std::cout << (*macroValuesTmp[level+1])(0)[204][0] << " " << (*macroValuesTmp[level+1])(0)[68][0] << std::endl;
      std::cout << (*macroR[level+1])(0)[204][0] << " " << (*macroR[level+1])(0)[68][0] << std::endl;
    }
#endif
    mySmoother->smooth(level+1, *macroValues[level+1], *macroValuesTmp[level+1],
                       numSmooths_pre[level+1]);
  }
  for(int level = num_levels; level > 0; --level) {
    multiGridLevels[level]->Prolong(*multiGridLevels[level-1], *macroValuesOld[level], *macroValues[level], *macroValues[level-1],prolong_relax_factor);
    mySmoother->smooth(level-1, *macroValues[level-1], *macroValuesTmp[level-1],
                       numSmooths_post[level-1]);
  }

  x = *macroValues[0];
}

template<class Scalar, int dim, class Scalar2>
void MultiGridKernel<Scalar,dim,Scalar2>::cycleW(int level,DistSVec<Scalar2,dim>& f, 
                                                 DistSVec<Scalar2,dim>& x) {

  mySmoother->smooth(level, x, f,numSmooths_pre[level]);

  if (level < num_levels) {

    multiGridLevels[level+1]->Restrict(*multiGridLevels[level], x,
                                       *macroValues[level+1]);

    // Restrict the residual. (I_h*d_h)
    multiGridLevels[level+1]->Restrict(*multiGridLevels[level], *macroR[level], *macroR[level+1]);
    *macroValuesOld[level+1] = *macroValues[level+1];
    if (isGeometric) {
      mySmoother->applyOperator(level+1, *macroValues[level+1], *macroValuesTmp[level+1]);
    } else {

      *macroValuesTmp[level+1] = 0.0;
    }
    
    *macroValuesTmp[level+1] += *macroR[level+1]*restrict_relax_factor;

    for (int mc = 0; mc < 2; ++mc) {

      cycleW(level+1, *macroValuesTmp[level+1], *macroValues[level+1]);
    }

    multiGridLevels[level+1]->Prolong(*multiGridLevels[level], *macroValuesOld[level+1], *macroValues[level+1], x,prolong_relax_factor);
  }

  mySmoother->smooth(level, x, f, numSmooths_post[level]);  
}

template<class Scalar, int dim, class Scalar2>
void MultiGridKernel<Scalar,dim,Scalar2>::smooth(int level,DistSVec<Scalar2,dim>& x,
                                               const DistSVec<Scalar2,dim>& f,int steps)
{
  
  static int l0_count = 0;
  if (level == 0)
    ++l0_count;
  static int r_tag = 0;
   /* 
  if (level == num_levels) {
    coarseSolver->setup(0,0,const_cast<DistSVec<Scalar2,dim>&>(f));
    coarseSolver->solve(const_cast<DistSVec<Scalar2,dim>&>(f),x);
    return;
  }
*/
  double rcurr;
  double scale = 1.0;
  MatVecProdH1<dim,Scalar2,dim>* mvp = dynamic_cast< MatVecProdH1<dim,Scalar2,dim>*>(macroA[level]);

  if (output) {
    domain->getCommunicator()->fprintf(stdout,"======================= Level %i ====================\n",level);
    domain->getCommunicator()->fprintf(stdout,"%d\n",l0_count);
  }
  for (int i = 0; i < steps; ++i) {
      
    double scale = 1.0;
    if (mvp) 
      mvp->apply(x, *macroR[level]);
    else
      multiGridLevels[level]->computeMatVecProd(*dynamic_cast< DistMvpMatrix<Scalar2,dim>*>(macroA[level]),
                                                x, *macroR[level]);

    (*macroR[level]) *= -1.0;
    (*macroR[level]) += f;

      if (level == 0) {
        DistSVec<double,1> tmp(macroR[level]->info());
        for (int l = 0; l < 1/*dim*/; ++l) {
          char fn[32];
#pragma omp parallel for
      for(int iSub = 0; iSub < numLocSub; ++iSub) {
          for (int k = 0; k < (*macroR[level])(iSub).size(); ++k)
            tmp(iSub)[k][0] = (*macroR[level])(iSub)[k][l];
      }
          sprintf(fn,"myR%d",l);
          domain->writeVectorToFile(fn,r_tag,(double)r_tag,tmp,
                                    &scale);
        }
        ++r_tag;
      }

      
    double rnew = (*macroR[level]).norm();
      // bailout?
      /*if (i > 0 && rnew > rcurr && x.norm() > 0.0) {
        x -= relaxationFactor*(*macroDX[level]);
        break;
      } else {*/
        rcurr = rnew;
//      }
     
    if (output)
      domain->getCommunicator()->fprintf(stdout,"i = %i r = %lf\n",i,(*macroR[level]).norm());

#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub) {
       smoothingMatrices[level][iSub]->smooth((*macroR[level])(iSub),(*macroDX[level])(iSub));
    }

    multiGridLevels[level]->assemble(*macroDX[level]);
    x += relaxationFactor*(*macroDX[level]);
  }
    
  // Recompute the residual
  if (mvp) 
    mvp->apply(x, *macroR[level]);
  else
    multiGridLevels[level]->computeMatVecProd(*dynamic_cast< DistMvpMatrix<Scalar2,dim>*>(macroA[level]),
                                              x, *macroR[level]);
  (*macroR[level]) *= -1.0;
  (*macroR[level]) += f;
      if (level == 0) {
        DistSVec<double,1> tmp(macroR[level]->info());
        for (int l = 0; l < 1/*dim*/; ++l) {
          char fn[32];
#pragma omp parallel for
      for(int iSub = 0; iSub < numLocSub; ++iSub) {
          for (int k = 0; k < (*macroR[level])(iSub).size(); ++k)
            tmp(iSub)[k][0] = (*macroR[level])(iSub)[k][l];
      }
          sprintf(fn,"myR%d",l);
          domain->writeVectorToFile(fn,r_tag,(double)r_tag,tmp,
                                    &scale);
        }
        ++r_tag;
      }
  if (output) {
    domain->getCommunicator()->fprintf(stdout,"final r = %lf\n",(*macroR[level]).norm());
  }
}

template<class Scalar, int dim, class Scalar2>
void MultiGridKernel<Scalar,dim,Scalar2>::applyOperator(int level,DistSVec<Scalar2,dim>& f,
                                                        DistSVec<Scalar2,dim>& x) {

  multiGridLevels[level]->computeMatVecProd(*dynamic_cast< DistMvpMatrix<Scalar2,dim>*>(macroA[level]),
                                            f,x);
}

template<class Scalar, int dim, class Scalar2>
void MultiGridKernel<Scalar,dim,Scalar2>::kspSolve(int lvl, DistSVec<Scalar2,dim>& f,
                                                   DistSVec<Scalar2,dim>& x) {

  coarseSolvers[lvl]->setup(0,0,const_cast<DistSVec<Scalar2,dim>&>(f));
  coarseSolvers[lvl]->solve(const_cast<DistSVec<Scalar2,dim>&>(f),x);
}

#define INSTANTIATION_HELPER(T,dim) \
    template class MultiGridKernel<T, dim, double>;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,5);
INSTANTIATION_HELPER(float,6);
INSTANTIATION_HELPER(float,7);

INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,5);
INSTANTIATION_HELPER(double,6);
INSTANTIATION_HELPER(double,7);
