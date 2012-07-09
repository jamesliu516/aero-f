#include <MultiGridPrec.h>

#include <Domain.h>
#include <DistGeoState.h>
#include <MultiGridLevel.h>
#include <MatVecProd.h>
#include <DistMvpMatrix.h>

#include <IoData.h>

//------------------------------------------------------------------------------

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
MultiGridPrec<Scalar,dim,Scalar2>::MultiGridPrec(Domain *dom, DistGeoState& distGeoState, PcData& pcData,KspData& coarseSolverData, bool createFineA, int **nodeType, BCApplier *bcs)
    : DistMat<Scalar2,dim>(dom), domain(dom), num_levels(pcData.num_multigrid_levels), agglom_size(8), numLocSub(dom->getNumLocSub()), multiGridLevels(new MultiGridLevel<Scalar2>*[num_levels+1]),
    macroValues(new DistSVec<Scalar2,dim>*[num_levels]), macroValuesTmp(new DistSVec<Scalar2,dim>*[num_levels]), geoState(distGeoState), pcData(pcData), initialized(false), distGeoState(distGeoState), coarseSolverData(coarseSolverData)
{

  ownsFineA = createFineA;

}

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::initialize() {

  initialized = true;

  smoothingMatrices = new MultiGridSmoothingMatrix<Scalar2,dim>**[num_levels+1];

  nSmooth1 = pcData.num_multigrid_smooth1;
  nSmooth2 = pcData.num_multigrid_smooth2;

  fine_sweeps = pcData.num_fine_sweeps;

  relaxationFactor = pcData.mg_smooth_relax;

  output = pcData.mg_output;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    domain->getSubDomain()[iSub]->getEdges().updateLength(distGeoState.getXn()(iSub));
    domain->getSubDomain()[iSub]->makeMasterFlag(domain->getNodeDistInfo());
  }

  multiGridLevels[0] = new MultiGridLevel<Scalar2>(*domain, domain->getNodeDistInfo(), domain->getEdgeDistInfo());
  multiGridLevels[0]->copyRefinedState(domain->getNodeDistInfo(), domain->getEdgeDistInfo(), geoState, *domain);

  char top_file_name[256];
  for(int level = 0; level < num_levels; ++level) {
    multiGridLevels[level+1] = new MultiGridLevel<Scalar2>(*domain, multiGridLevels[level]->getNodeDistInfo(), multiGridLevels[level]->getEdgeDistInfo());
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
                                          NULL);
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
  
  for(int level = 0; level <= num_levels; ++level) {
    macroA[level] = 0;
    if (level > 0 || ownsFineA)
      macroA[level] = new DistMvpMatrix<Scalar2,dim>(domain,multiGridLevels[level]->getNodeDistInfo().subLen,
                                                  multiGridLevels[level]->getEdges(),
                                                  NULL); 

    macroValues[level] = new DistSVec<Scalar2,dim>(multiGridLevels[level]->getNodeDistInfo());
    macroValuesOld[level] = new DistSVec<Scalar2,dim>(multiGridLevels[level]->getNodeDistInfo());
    macroValuesTmp[level] = new DistSVec<Scalar2,dim>(multiGridLevels[level]->getNodeDistInfo());
    macroR[level] = new DistSVec<Scalar2,dim>(multiGridLevels[level]->getNodeDistInfo());
    macroIrey[level] = new DistVec<Scalar2>(multiGridLevels[level]->getNodeDistInfo());
    macroDX[level] = new DistSVec<Scalar2,dim>(multiGridLevels[level]->getNodeDistInfo());
  }

  typename MultiGridSmoothingMatrix<Scalar2,dim>::SmoothingMode smoothing_mode;
  if (pcData.mg_smoother == PcData::MGJACOBI)
    smoothing_mode = MultiGridSmoothingMatrix<Scalar2,dim>::BlockJacobi;
  else if (pcData.mg_smoother == PcData::MGLINEJACOBI)
    smoothing_mode = MultiGridSmoothingMatrix<Scalar2,dim>::LineJacobi;
  else if (pcData.mg_smoother == PcData::MGRAS)
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

  coarseMvp = new MultiGridPrecMatVecProd<Scalar2,dim>(macroA[num_levels], 
                                                      multiGridLevels[num_levels]);

  coarsePrec = new MultiGridPrecJacobiPrec<Scalar2,dim>(smoothingMatrices[num_levels],
                                                        multiGridLevels[num_levels]);

  coarseSolver = new GmresSolver<DistSVec<Scalar2,dim>, MultiGridPrecMatVecProd<Scalar2,dim>,
                                 MultiGridPrecJacobiPrec<Scalar2,dim>, Communicator>(
                     macroValues[num_levels]->info(), coarseSolverData, coarseMvp,
                     coarsePrec, domain->getCommunicator());

  coarseSolver->setEps(1.0e-12);
  coarseSolver->setMaxIts(100);
  coarseSolver->disableOutput();

    double scale = 1.0;
    DistSVec<double,1> tmp(macroR[0]->info());
#pragma omp parallel for
      for(int iSub = 0; iSub < numLocSub; ++iSub) {

        for (int k = 0; k < (*macroR[0])(iSub).size(); ++k)
          tmp(iSub)[k][0] = static_cast<double>(/*domain->getSubDomain()[iSub]->getNodeMap()[*/(multiGridLevels[1]->getNodeMapping())(iSub)[k]/*]*/);
     }   
      char fn[32];
     
   sprintf(fn,"nodeMapping");
        domain->writeVectorToFile(fn,0,(double)0.0,tmp,
                                  &scale);
      
  
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
MultiGridPrec<Scalar,dim,Scalar2>::~MultiGridPrec()
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

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::setup()
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

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::apply(DistSVec<Scalar2,dim> & x, DistSVec<Scalar2,dim> & Px)
{
  *macroValues[0] = 0.0;
  *macroValues[1] = 0.0;
  *macroValuesOld[1] = 0.0;//*macroValues[level+1];
  //multiGridLevels[1]->Prolong(*multiGridLevels[0], *macroValuesOld[1], *macroValues[1], x);
  smooth(0, *macroValues[0], x,nSmooth1);
  *macroValuesTmp[0] = x;
  *macroValuesOld[0] = 0.0;//*macroValues[0];
   MatVecProdH1<dim,Scalar2,dim>* mvp = dynamic_cast< MatVecProdH1<dim,Scalar2,dim>*>(macroA[0]);
//  mvp->apply(*macroValues[0], Px);
//  std::cout << "fn 0 " << Px*Px << std::endl;
  /*char filename[256];
  sprintf(filename,"result_0_1.xpost");
  multiGridLevels[0]->writeXpostFile(filename, *macroValues[0],0);
  sprintf(filename,"R_0_1.xpost");
  multiGridLevels[0]->writeXpostFile(filename, *macroR[0],0);
  exit(-1);
 */
  for(int level = 0; level < num_levels; ++level) {
    *macroValuesOld[level+1] = 0.0;//*macroValues[level+1];
    *macroValues[level+1] = 0.0;
    //multiGridLevels[level+1]->Prolong(*multiGridLevels[level], *macroValuesOld[level+1], *macroValues[level+1], *macroValues[level]);
    //multiGridLevels[level+1]->Restrict(*multiGridLevels[level], *macroValues[level], *macroValues[level+1]);
    // Restrict the residual. (I_h*d_h)
    multiGridLevels[level+1]->Restrict(*multiGridLevels[level], *macroR[level], *macroR[level+1]);
    //multiGridLevels[level]->assembleMax(*macroR[level]);
    //*macroValues[level+1] = 0.0;
    *macroValuesOld[level+1] = 0.0;//*macroValues[level+1];
    // Now compute A_H*u_H
/*    multiGridLevels[level+1]->computeMatVecProd(*dynamic_cast< DistMvpMatrix<Scalar2,dim>*>(macroA[level+1]),
                                                *macroValues[level+1], *macroValuesTmp[level+1]);
    */
//    std::cout << "fn " << level+1 << " " << pow(macroValuesTmp[level+1]->norm(),2.0) << std::endl;
    *macroValuesTmp[level+1] = *macroR[level+1];
    //*macroValues[level+1] = 0.01; 
    smooth(level+1, *macroValues[level+1], *macroValuesTmp[level+1],
           ((level+1==num_levels)?fine_sweeps:nSmooth1));
    //multiGridLevels[level+1]->assembleMax(*macroValues[level+1]);
    //multiGridLevels[level+1]->assembleMax(*macroValues[level+1]);
    /*sprintf(filename,"macroR.xpost",level+1);
    multiGridLevels[0]->writeXpostFile(filename, *macroValues[level+1],0);
    exit(-1);*/
    /*sprintf(filename,"result_%i_1.xpost",level+1);
    multiGridLevels[0]->writeXpostFile(filename, *macroValues[level+1],0);
    sprintf(filename,"R_%i_1.xpost",level+1);
    multiGridLevels[0]->writeXpostFile(filename, *macroValuesTmp[level+1],0);
    */
  }
  for(int level = num_levels; level > 0; --level) {
    //multiGridLevels[level]->assembleMax(*macroValues[level]);
    //multiGridLevels[0]->assembleMax(*macroValues[0]);
    multiGridLevels[level]->Prolong(*multiGridLevels[level-1], *macroValuesOld[level], *macroValues[level], *macroValues[level-1]);
    //multiGridLevels[level-1]->assembleMax(*macroValues[level-1]);
    //sprintf(filename,"prolonged.xpost");
    //multiGridLevels[0]->writeXpostFile(filename, *macroValues[level-1],0);
    //*macroValuesTmp[level-1] = *macroValues[level-1];
    //multiGridLevels[0]->assembleMax(*macroValues[0]);
    smooth(level-1, *macroValues[level-1], *macroValuesTmp[level-1],nSmooth2);
  }
  //exit(-1);
  /*sprintf(filename,"final.xpost");
  multiGridLevels[0]->writeXpostFile(filename, *macroValues[0],0);
  exit(-1);*/
  //smooth(0, *macroValues[0], x);
  //multiGridLevels[0]->assembleMax(*macroValues[0]);
  Px = *macroValues[0];
}

//------------------------------------------------------------------------------
template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::smooth(int level,DistSVec<Scalar2,dim>& x,
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
  MatVecProdH1<dim,Scalar2,dim>* mvp = dynamic_cast< MatVecProdH1<dim,Scalar2,dim>*>(macroA[level]);

  if (output) {
    domain->getCommunicator()->fprintf(stdout,"======================= Level %i ====================\n",level);
    domain->getCommunicator()->fprintf(stdout,"%d\n",l0_count);
  }
  //MatVecProdH1<dim,Scalar2,dim>* mvp = dynamic_cast< MatVecProdH1<dim,Scalar2,dim>*>(macroA[level]);
  /*if(mvp == NULL) {
      // std::cout << "Error: multigrid preconditioner must be used with an H1 product" << std::endl;
      x = f;
  } else {*/
    for (int i = 0; i < steps; ++i) {
      // std::cout << "Smooth step " << i << std::endl;

      // Compute the residual
      
      double scale = 1.0;
      if (mvp) 
        mvp->apply(x, *macroR[level]);
      else
        multiGridLevels[level]->computeMatVecProd(*dynamic_cast< DistMvpMatrix<Scalar2,dim>*>(macroA[level]),
                                                  x, *macroR[level]);

      (*macroR[level]) *= -1.0;
      (*macroR[level]) += f;
      
      if (level == 0 && l0_count >= 0) {
        DistSVec<double,1> tmp(macroR[level]->info());
        for (int l = 0; l < dim; ++l) {
#pragma omp parallel for
          for(int iSub = 0; iSub < numLocSub; ++iSub) {

            for (int k = 0; k < (*macroR[level])(iSub).size(); ++k)
              tmp(iSub)[k][0] = (*macroR[level])(iSub)[k][l];
          }
          char fn[32];
          sprintf(fn,"myR%d",l);
          domain->writeVectorToFile(fn,r_tag,(double)r_tag,tmp,
                                    &scale);
#pragma omp parallel for
      for(int iSub = 0; iSub < numLocSub; ++iSub) {
          for (int k = 0; k < (*macroR[level])(iSub).size(); ++k)
            tmp(iSub)[k][0] = (x)(iSub)[k][l];
      }         
          sprintf(fn,"myU%d",l);
          domain->writeVectorToFile(fn,r_tag,(double)r_tag,tmp,
                                    &scale);
        }
        ++r_tag;
      }

/*
      if (l0_count == 2 && i == 0)
        multiGridLevels[level]->writeXpostFile("RafterProlong.xpost",*macroR[level],0);
  */      

      double rnew = (*macroR[level]).norm();
      // bailout?
      /*if (i > 0 && rnew > rcurr && x.norm() > 0.0) {
        x -= relaxationFactor*(*macroDX[level]);
        break;
      } else {*/
        rcurr = rnew;
//      }
     
      //DistSVec<Scalar2,dim> tmp(f); 
      //multiGridLevels[level]->assembleMax(tmp);//*macroR[level]);

      if (output)
        domain->getCommunicator()->fprintf(stdout,"i = %i r = %lf\n",i,(*macroR[level]).norm());
     // std::cout << "dot = " << (*macroR[level])*f / (sqrt((*macroR[level]).norm() * f.norm())) << std::endl;

#pragma omp parallel for
      for(int iSub = 0; iSub < numLocSub; ++iSub) {
        smoothingMatrices[level][iSub]->smooth((*macroR[level])(iSub),(*macroDX[level])(iSub));
      }

      multiGridLevels[level]->assemble(*macroDX[level]);
      // std::cout << "DX = " << (*macroDX[level]).norm() << std::endl;

      // *macroDX[level] = *macroR[level];
      x += relaxationFactor*(*macroDX[level]);
      /*if (level == 0) 
        ++r_tag;
      if (r_tag == 40)
        exit(-1);
      */
    }
      

    // Recompute the residual
    if (mvp) 
      mvp->apply(x, *macroR[level]);
    else
      multiGridLevels[level]->computeMatVecProd(*dynamic_cast< DistMvpMatrix<Scalar2,dim>*>(macroA[level]),
                                                x, *macroR[level]);
    (*macroR[level]) *= -1.0;
    (*macroR[level]) += f;
    if (output) {
      domain->getCommunicator()->fprintf(stdout,"final r = %lf\n",(*macroR[level]).norm());
    }
    
  if (level == 0 && l0_count >= 0) {
    double scale = 1.0;
    DistSVec<double,1> tmp(macroR[level]->info());
    for (int l = 0; l < dim; ++l) {
#pragma omp parallel for
      for(int iSub = 0; iSub < numLocSub; ++iSub) {

        for (int k = 0; k < (*macroR[level])(iSub).size(); ++k)
          tmp(iSub)[k][0] = (*macroR[level])(iSub)[k][l];
      } 
        char fn[32];
        sprintf(fn,"myR%d",l);
        domain->writeVectorToFile(fn,r_tag,(double)r_tag,tmp,
                                    &scale);
#pragma omp parallel for
      for(int iSub = 0; iSub < numLocSub; ++iSub) {
        for (int k = 0; k < (*macroR[level])(iSub).size(); ++k)
          tmp(iSub)[k][0] = (x)(iSub)[k][l];
      }  
        sprintf(fn,"myU%d",l);
        domain->writeVectorToFile(fn,r_tag,(double)r_tag,tmp,
                                  &scale);
      }
      ++r_tag;
    }
  
  //}
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::getData(DistMat<Scalar2,dim>& mat)/* ,
                                                DistSVec<Scalar2,dim>& V,
                                                SpaceOperator<dim>& spo,
                                                DistTimeState<dim>* timeState)*/
{
//  VarFcn* varFcn = spo.getVarFcn();

  macroA[0] = &mat;
  
/*  *macroValues[0] = V;
  if (timeState)
    *macroIrey[0] = *timeState->getInvReynolds();
  else
    *macroIrey[0] = 0.0;
*/
  //std::cout << "Vnorm 0 = " << V*V <<std::endl;
  //smooth(0, *macroValues[0], V);
/*#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    for (int i = 0; i < macroValues[0]->subSize(iSub); ++i)
      std::cout << i << " " << (*macroValues[0])(iSub)[i][0] << " " << 
        (*timeState)(iSub).getDt()[i] << std::endl;
  }
*/
}

//------------------------------------------------------------------------------
#define INSTANTIATION_HELPER(T,dim) \
    template class MultiGridPrec<T, dim, double>;

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
