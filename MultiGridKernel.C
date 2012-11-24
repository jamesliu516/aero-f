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

template<class Scalar>
MultiGridKernel<Scalar>::MultiGridKernel(Domain *dom, DistGeoState& distGeoState, IoData& ioData,int num_levels)
    :  domain(dom), num_levels(num_levels), agglom_size(8), numLocSub(dom->getNumLocSub()), multiGridLevels(new MultiGridLevel<Scalar>*[num_levels]),
    initialized(false),ioData(ioData), geoState(distGeoState)
{

  beta = ioData.mg.directional_coarsening_factor;

  fixLocations = new std::set<int>*[num_levels];
  for (int lvl = 0; lvl < num_levels; ++lvl) {

    fixLocations[lvl] = new std::set<int>[numLocSub];
  }

  coarsen4to1 = (ioData.mg.coarseningRatio == MultiGridData::FOURTOONE);
}

template<class Scalar>
void MultiGridKernel<Scalar>::setParameters(int v1, int v2, int
finesweeps, double relax, int do_out) {

  nSmooth1 = v1;
  nSmooth2 = v2;
  relaxationFactor = relax;
  fine_sweeps = finesweeps;
  output = do_out;
}

template<class Scalar>
void MultiGridKernel<Scalar>::initialize(int dim,int neq1,int neq2) {

  initialized = true;

  MultiGridMethod mgm = MultiGridGeometric;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    domain->getSubDomain()[iSub]->getEdges().updateLength(geoState.getXn()(iSub));
    domain->getSubDomain()[iSub]->makeMasterFlag(domain->getNodeDistInfo());
  }

  multiGridLevels[0] = new MultiGridLevel<Scalar>(mgm,NULL,*domain, domain->getNodeDistInfo(), domain->getEdgeDistInfo());
  multiGridLevels[0]->copyRefinedState(domain->getNodeDistInfo(), domain->getEdgeDistInfo(),domain->getInletNodeDistInfo(), domain->getFaceDistInfo(),geoState, *domain,dim,neq1,neq2);

  int maxNumNodesPerAgglom = 6;
  if (coarsen4to1)
    maxNumNodesPerAgglom = 3;

  char top_file_name[256];
  for(int level = 0; level < num_levels-1; ++level) {
    multiGridLevels[level+1] = new MultiGridLevel<Scalar>(mgm,multiGridLevels[level],*domain, multiGridLevels[level]->getNodeDistInfo(), multiGridLevels[level]->getEdgeDistInfo());
    multiGridLevels[level+1]->agglomerate(multiGridLevels[level]->getNodeDistInfo(),
                                          multiGridLevels[level]->getEdgeDistInfo(),
                                          multiGridLevels[level]->getIdPat(),
                                          multiGridLevels[level]->getSharedNodes(),
                                          multiGridLevels[level]->getConnectivity(),
                                          multiGridLevels[level]->getEdges(),
                                          multiGridLevels[level]->getSharedEdges(),
                                          multiGridLevels[level]->getNumSharedEdges(),
                                          *domain,dim,neq1,neq2,
                                          multiGridLevels[level]->getEdgeNormals(),
                                          multiGridLevels[level]->getCtrlVol(),
                                          NULL,beta, maxNumNodesPerAgglom);

    if (coarsen4to1) {

      MultiGridLevel<Scalar>* mg2 = new MultiGridLevel<Scalar>(mgm,multiGridLevels[level+1],*domain, multiGridLevels[level+1]->getNodeDistInfo(), multiGridLevels[level+1]->getEdgeDistInfo());
      mg2->agglomerate(multiGridLevels[level+1]->getNodeDistInfo(),
                       multiGridLevels[level+1]->getEdgeDistInfo(),
                       multiGridLevels[level+1]->getIdPat(),
                       multiGridLevels[level+1]->getSharedNodes(),
                       multiGridLevels[level+1]->getConnectivity(),
                       multiGridLevels[level+1]->getEdges(),
                       multiGridLevels[level+1]->getSharedEdges(),
                       multiGridLevels[level+1]->getNumSharedEdges(),
                       *domain,dim,neq1,neq2,
                       multiGridLevels[level+1]->getEdgeNormals(),
                       multiGridLevels[level+1]->getCtrlVol(),
                       NULL,beta, maxNumNodesPerAgglom);
      mg2->mergeFinerInto(*multiGridLevels[level+1]);
      //delete multiGridLevels[level+1];
      multiGridLevels[level+1] = mg2;
    }

    sprintf(top_file_name,"level%d",level+1);
    multiGridLevels[level+1]->writePVTUFile(top_file_name);
    sprintf(top_file_name,"agglevel%d",level+1);
    multiGridLevels[level+1]->writePVTUAgglomerationFile(top_file_name);
    domain->getCommunicator()->fprintf(stdout,"Agglomerated level %d\n", level+1);
    fflush(stdout);
  }

}

template<class Scalar>
void MultiGridKernel<Scalar>::setUseVolumeWeightedAverage(bool b) {

  for(int level = 0; level < num_levels; ++level) {

    multiGridLevels[level]->setUseVolumeWeightedAverage(b);
  }
}

template<class Scalar>
MultiGridKernel<Scalar>::~MultiGridKernel()
{
  for (int level = 0; level < num_levels; ++level) {
    delete multiGridLevels[level];
  }
  delete []multiGridLevels;
}

//------------------------------------------------------------------------------

template <class Scalar>
template<class Scalar2, int dim>
void MultiGridKernel<Scalar>::Restrict(int coarseLvl, DistSVec<Scalar2,dim>& fine,
                                       DistSVec<Scalar2,dim>& coarse) {

  multiGridLevels[coarseLvl]->Restrict(*multiGridLevels[coarseLvl-1],
                                       fine,coarse);
}

template <class Scalar>
template<class Scalar2, int dim>
void MultiGridKernel<Scalar>::Prolong(int coarseLvl,
                                      DistSVec<Scalar2,dim>& coarseOld,
                                      DistSVec<Scalar2,dim>& coarse,
                                      DistSVec<Scalar2,dim>& fine,double relax) {

  multiGridLevels[coarseLvl]->Prolong(*multiGridLevels[coarseLvl-1],
                                      coarseOld,coarse,fine,relax);
}

template <class Scalar>
template<class Scalar2, int dim>
void MultiGridKernel<Scalar>::
fixNegativeValues(int lvl,DistSVec<Scalar2,dim>& V, 
                  DistSVec<Scalar2,dim>& U, 
                  DistSVec<Scalar2,dim>& dx, 
                  DistSVec<Scalar2,dim>& f, 
                  DistSVec<Scalar2,dim>& forig, 
                  VarFcn* vf,
                  MultiGridOperator<Scalar2,dim>* op) {

  int rnk;
  MPI_Comm_rank(MPI_COMM_WORLD,&rnk);

#pragma omp parallel for 
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
  
    SVec<Scalar2,dim>& Vl = V(iSub);
    SVec<Scalar2,dim>& Ul = U(iSub);
    SVec<Scalar2,dim>& dxl = dx(iSub);
    SVec<Scalar2,dim>& fl = f(iSub);
    SVec<double,3>& X = multiGridLevels[lvl]->getXn()(iSub);
    DistGeoState& geoState = multiGridLevels[lvl]->getGeoState();

    DistVec<Vec3D>& edgeNorm = geoState.getEdgeNormal();
    DistVec<double>& ctrlVol = geoState.getCtrlVol();

    for (int i = 0; i < Vl.size(); ++i) {

      if (vf->getPressure(Vl[i]) <= 0.0 ||
          vf->getDensity(Vl[i]) <= 0.0 /*|| (iSub == 0 && rnk == 17 && i == 1604 && lvl == 1)*/) {

        fprintf(stderr,"lvl = %d iSub = %d, rank = %d\n", lvl,iSub, rnk);
        fprintf(stderr,"Fixed negative value (p, rho) = (%lf, %lf) at"
                       " node %d (x = [%lf, %lf, %lf])\n", vf->getPressure(Vl[i]), vf->getDensity(Vl[i]),
                       i, X[i][0],X[i][1],X[i][2]);

        Connectivity* con = multiGridLevels[lvl]->getConnectivity()[iSub];
        fprintf(stderr,"f[i] = [%e, %e, %e, %e, %e]\n",
                fl[i][0], fl[i][1], fl[i][2], fl[i][3], fl[i][4]);
        fprintf(stderr,"U[i] = [%e, %e, %e, %e, %e]\n",
                Ul[i][0]-dxl[i][0], Ul[i][1]-dxl[i][1], Ul[i][2]-dxl[i][2], Ul[i][3]-dxl[i][3], Ul[i][4]-dxl[i][4]);
        fprintf(stderr,"dt = %lf\n", op->queryTimeStep(iSub,i));
        fprintf(stderr,"nodetype[i] = %d\n", multiGridLevels[lvl]->getNodeType(iSub)[i]);
        fprintf(stderr,"ctrlvol[i] = %e\n", ctrlVol(iSub)[i]);
        for (int j = 0; j < con->num(i); ++j) {

          if ((*con)[i][j] == i) continue;
          int l = (*con)[i][j];
          fprintf(stderr,"f[%d] = [%e, %e, %e, %e, %e]\n",l,
                  fl[l][0], fl[l][1], fl[l][2], fl[l][3], fl[l][4]);
          fprintf(stderr,"U[%d] = [%e, %e, %e, %e, %e]\n",l,
                Ul[l][0]-dxl[l][0], Ul[l][1]-dxl[l][1], Ul[l][2]-dxl[l][2], Ul[l][3]-dxl[l][3], Ul[l][4]-dxl[l][4]);
          fprintf(stderr,"dt = %lf\n", op->queryTimeStep(iSub,l));
          fprintf(stderr,"nodetype[%d] = %d\n",l, multiGridLevels[lvl]->getNodeType(iSub)[l]);
          int edge_l = multiGridLevels[lvl]->getEdges()[iSub]->findOnly(i,l);
          fprintf(stderr,"Edge_normal = [%e %e %e]; area = %e\n", edgeNorm(iSub)[edge_l][0],edgeNorm(iSub)[edge_l][1],edgeNorm(iSub)[edge_l][2],edgeNorm(iSub)[edge_l].norm());
        fprintf(stderr,"ctrlvol[%d] = %e\n", l, ctrlVol(iSub)[l]);


        }
 
        if (vf->getPressure(Vl[i]) > 0.0 && vf->getDensity(Vl[i]) > 0.0 ) continue;

        for (int k = 0; k < dim; ++k) {
          Ul[i][k] -= dxl[i][k];
          fl[i][k] = forig[i][k];
        }
        vf->conservativeToPrimitive(Ul[i],Vl[i]);
        fixLocations[lvl][iSub].insert(i);
      }
    }
  }
}

template <class Scalar>
template<class Scalar2, int dim>
void MultiGridKernel<Scalar>::
applyFixes(int lvl,DistSVec<Scalar2,dim>& f) {

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    for (std::set<int>::const_iterator it = fixLocations[lvl][iSub].begin();
         it != fixLocations[lvl][iSub].end(); ++it) {
      for (int k = 0; k < dim; ++k) {
        f(iSub)[*it][k] = 0.0;
      }
    }
  }
}
 

#define INSTANTIATION_HELPER(T) \
    template class MultiGridKernel<T>;

#define INSTANTIATION_HELPER2(T,T2,D) \
 template void MultiGridKernel<T>::Restrict(int coarseLvl, DistSVec<T2,D>&, \
                                   DistSVec<T2,D>&); \
 template void MultiGridKernel<T>::Prolong(int coarseLvl, DistSVec<T2,D>&, \
                                   DistSVec<T2,D>&,DistSVec<T2,D>&, \
                                   double); \
template void MultiGridKernel<T>:: \
fixNegativeValues(int,DistSVec<T2,D>& V, \
                  DistSVec<T2,D>& U,  \
                  DistSVec<T2,D>& dx,  \
                  DistSVec<T2,D>& f, \
                  DistSVec<T2,D>& forig, \
                  VarFcn* vf, MultiGridOperator<T2,D>*);\
template void MultiGridKernel<T>:: \
applyFixes(int,DistSVec<T2,D>& V);



//INSTANTIATION_HELPER(float);

INSTANTIATION_HELPER(double);
INSTANTIATION_HELPER2(double,double,1);
INSTANTIATION_HELPER2(double,double,2);
INSTANTIATION_HELPER2(double,double,5);
INSTANTIATION_HELPER2(double,double,6);
INSTANTIATION_HELPER2(double,double,7);
