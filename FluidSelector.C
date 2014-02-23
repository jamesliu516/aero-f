#include "LevelSet/LevelSetStructure.h"
//#include <Domain.h>
#include <cassert>

//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::initializeFluidIds(DistSVec<double,dim> &Phin, DistSVec<double,dim> &Phinm1, DistSVec<double,dim> &Phinm2){
  getFluidId(Phin);
  *fluidIdn = *fluidId;
  if(fluidIdnm1) getFluidId(*fluidIdnm1, Phinm1);
  if(fluidIdnm2) getFluidId(*fluidIdnm2, Phinm2);
}

//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::reinitializeFluidIds(DistVec<int> &fsId, DistSVec<double,dim> &Phin)
{
  getFluidId(*fluidId, Phin, &fsId);
  //if (fluidIdnm1) getFluidId(*fluidIdm1, Phinm1, &fsId);
  //if(fluidIdnm2) getFluidId(*fluidIdm2, Phinm1, &fsId);
  *fluidIdn = *fluidId;
  if(fluidIdnm1) *fluidIdnm1 = *fluidId;
  if(fluidIdnm2) *fluidIdnm2 = *fluidId;
}

template<int dim>
void FluidSelector::reinitializeFluidIdsWithCracking(DistVec<int> &fsId, DistSVec<double,dim> &Phin)
{
  getFluidId(*fluidId, Phin, NULL);
  //if (fluidIdnm1) getFluidId(*fluidIdm1, Phinm1, &fsId);
  //if(fluidIdnm2) getFluidId(*fluidIdm2, Phinm1, &fsId);
  *fluidIdn = *fluidId;
  if(fluidIdnm1) *fluidIdnm1 = *fluidId;
  if(fluidIdnm2) *fluidIdnm2 = *fluidId;
}


//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::getFluidId(DistSVec<double,dim> &Phi){
  assert(dim<=numPhases-1);
  int numLocSub = Phi.numLocSub();
  int iSub;
#pragma omp parallel for
  for(iSub=0; iSub<numLocSub; ++iSub) {
    double (*phi)[dim] = Phi.subData(iSub);
    int     *tag       = fluidId->subData(iSub);
    int burnTag;
    for(int iNode=0; iNode<Phi.subSize(iSub); iNode++){
      //if (programmedBurn && programmedBurn->isBurnedEOS(tag[iNode],burnTag))
      //	continue;
      tag[iNode] = 0;
      for(int i=0; i<dim; i++) {
        if(phi[iNode][i]>0.0) {
	  if (programmedBurn && (programmedBurn->isUnburnedEOS(i+1,burnTag) ||
				 programmedBurn->isBurnedEOS(i+1,burnTag)) ) {
	    if (programmedBurn->nodeInside(burnTag,iSub,iNode) || 
                programmedBurn->isFinished(burnTag))
	      tag[iNode] = programmedBurn->getBurnedEOS(burnTag);
	    else
	      tag[iNode] = programmedBurn->getUnburnedEOS(burnTag);
	    break;
	  }
	  else 
	    {  tag[iNode] = i+1; break; }
	}
      }
    }
  }
}
//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::getFluidId(int &tag, double *phi){
  assert(dim<=numPhases-1);
  tag = 0;
  for(int i=0; i<dim; i++)
    if(phi[i]>0.0) {tag = i+1; return; }
}

//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::getFluidId(Vec<int> &tag, SVec<double,dim> &phi){
  assert(dim<=numPhases-1);
  int burnTag;
  for(int iNode=0; iNode<phi.size(); iNode++){
    tag[iNode] = 0;
    for(int i=0; i<dim; i++) {
      if(phi[iNode][i]>0.0) {
	if (programmedBurn && (programmedBurn->isUnburnedEOS(i+1,burnTag) ||
			       programmedBurn->isBurnedEOS(i+1,burnTag)) ) {
	  if (programmedBurn->nodeInside(burnTag,iNode) ||
              programmedBurn->isFinished(burnTag))
	    tag[iNode] = programmedBurn->getBurnedEOS(burnTag);
	  else
	    tag[iNode] = programmedBurn->getUnburnedEOS(burnTag);
	  break;
	}
	else
	  { tag[iNode] = i+1; break; }
      }
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::getFluidId(DistVec<int> &Tag, DistSVec<double,dim> &Phi, DistVec<int>* fsId){
  assert(dim<=numPhases-1);
  //std::cout << "Dim = " << dim << std::endl;
  int numLocSub = Phi.numLocSub();
  int oldtag;
  int iSub;
  //std::cout << programmedBurn << std::endl;
#pragma omp parallel for
  for(iSub=0; iSub<numLocSub; ++iSub) {
    double (*phi)[dim] = Phi.subData(iSub);
    int     *tag       = Tag.subData(iSub);
    int     *fsid      = fsId ? fsId->subData(iSub) : 0;
    int burnTag;
    for(int iNode=0; iNode<Phi.subSize(iSub); iNode++){
      //if (programmedBurn && programmedBurn->isBurnedEOS(tag[iNode],burnTag))
      //	continue;

      tag[iNode] = 0;
      if(fsid && fsid[iNode]!=0) //isolated by structure. use fsId as fluidId.
        tag[iNode] = fsid[iNode];
      else for(int i=0; i<dim; i++) {
	if(phi[iNode][i]>0.0) { 
	  if (programmedBurn && (programmedBurn->isUnburnedEOS(i+1,burnTag) ||
				 programmedBurn->isBurnedEOS(i+1,burnTag)) ) {
	    if (programmedBurn->nodeInside(burnTag,iSub,iNode) ||
                programmedBurn->isFinished(burnTag))
	      tag[iNode] = programmedBurn->getBurnedEOS(burnTag);
	    else
	      tag[iNode] = programmedBurn->getUnburnedEOS(burnTag);
	    break;
	  }
	  else{
	    tag[iNode] = i+1; break; 
	  }
	}
      }
    }
  }
}

//------------------------------------------------------------------------------

template<int dim> /*this dim is actually dimLS*/
void FluidSelector::updateFluidIdFS(DistLevelSetStructure *distLSS, DistSVec<double,dim> &PhiV)
{
  assert(dim<=numPhases-1);
  DistVec<int> &fsId(distLSS->getStatus());

  int burnTag;
  int iSub;
#pragma omp parallel for
  for (iSub=0; iSub<PhiV.numLocSub(); ++iSub) {
    Vec<int> &subfsId(fsId(iSub));
    Vec<int> &subId((*fluidId)(iSub));
    SVec<double,dim> &subPhiV(PhiV(iSub));
    LevelSetStructure &LSS((*distLSS)(iSub));

    for(int i=0; i<subPhiV.size(); i++) {
      int Id = LSS.fluidModel(0.0,i);
      bool swept = LSS.isSwept(0.0,i);
      if(subId[i]==2) fprintf(stderr,"my subId = %d, swept = %d, Id = %d.\n", subId[i], swept, Id);

      if(Id==0) { // not isolated by structure. need to consider level-set
        if(swept) {
          subId[i] = 0;
          for(int k=0; k<dim; k++) {
            if(subPhiV[i][k]>0.0) {
	      if (programmedBurn && (programmedBurn->isUnburnedEOS(k+1,burnTag) ||
				     programmedBurn->isBurnedEOS(k+1,burnTag)) ) {
		if (programmedBurn->nodeInside(burnTag,iSub,i) ||
                    programmedBurn->isFinished(burnTag)) {
                  fprintf(stderr,"Inside updateFluidIdFS, I am here!\n");
		  subId[i] = programmedBurn->getBurnedEOS(burnTag);}
		else
		  subId[i] = programmedBurn->getUnburnedEOS(burnTag);
		break;
	      }
	      else {
		subId[i] = k+1;
		break;
	      }
            }
	  }
	} 
      } else // isolated by structure. Id determined by intersector
        subId[i] = Id;
    }


    for(int i=0; i<subPhiV.size(); i++) {
      if(subId[i]==2) fprintf(stderr,"In FluidSelector::updateFluidIdFS(...), found Id = 2!\n");}
  }
}

//------------------------------------------------------------------------------

template<int dim> /*this dim is actually dimLS*/
void FluidSelector::updateFluidIdFS2(DistLevelSetStructure *distLSS, DistSVec<double,dim> &PhiV)
{
  if(programmedBurn) {fprintf(stderr,"ERROR: function 'updateFluidIdFS2' does not support Programmed Burn at the moment!\n");exit(-1);}
  //domain->updateFluidIdFS2(*distLSS, PhiV, *fluidId);
}

//------------------------------------------------------------------------------

template<int dim> /*this dim is actually dimLS*/
void FluidSelector::updateFluidIdFF(DistLevelSetStructure *distLSS, DistSVec<double,dim> &Phi)
{
  assert(dim<=numPhases-1);
  int numLocSub = Phi.numLocSub();
  int burnTag;
  DistVec<int> &fsId(distLSS->getStatus());
  int iSub;
#pragma omp parallel for
  for(iSub=0; iSub<numLocSub; ++iSub) {
    double (*phi)[dim] = Phi.subData(iSub);
    int     *tag       = fluidId->subData(iSub);
    int     *fsid      = fsId.subData(iSub);
    for(int iNode=0; iNode<Phi.subSize(iSub); iNode++){
      if(fsid[iNode]!=0) {
        if(fsid[iNode]!=tag[iNode]){
          fprintf(stderr,"This must be a bug!\n"); exit(-1);}
        continue;
      }
      tag[iNode] = 0;
      for(int i=0; i<dim; i++) {
	if(phi[iNode][i]>0.0) {
	  if (programmedBurn && (programmedBurn->isUnburnedEOS(i+1,burnTag) ||
				 programmedBurn->isBurnedEOS(i+1,burnTag)) ) {
	    if (programmedBurn->nodeInside(burnTag,iSub,iNode) ||
                programmedBurn->isFinished(burnTag)) {
	      tag[iNode] = programmedBurn->getBurnedEOS(burnTag);
              fprintf(stderr,"I am here... \n");
            }
	    else {
	      tag[iNode] = programmedBurn->getUnburnedEOS(burnTag);
            }
	    break;
	  }
	  else{ 
	    tag[iNode] = i+1; break; 
	  }
	}
      }
    }

    for(int iNode=0; iNode<Phi.subSize(iSub); iNode++)
      if(tag[iNode] == 2) fprintf(stderr," Caught Id = 2.\n");
  }
}

//------------------------------------------------------------------------------
// This function allows fracture. It does not look at fsid.
template<int dim> /*this dim is actually dimLS*/
void FluidSelector::updateFluidIdFF2(DistLevelSetStructure *distLSS, DistSVec<double,dim> &Phi)
{
  int numLocSub = Phi.numLocSub();
  int burnTag;
  int iSub;
#pragma omp parallel for
  for(iSub=0; iSub<numLocSub; ++iSub) {
    double (*phi)[dim]     = Phi.subData(iSub);
    int     *tag           = fluidId->subData(iSub);
    LevelSetStructure &LSS = (*distLSS)(iSub);

    for(int iNode=0; iNode<Phi.subSize(iSub); iNode++){
      if(LSS.isOccluded(0.0,iNode)) {
        //phi[iNode][dim-1] = 0.0;
        tag[iNode] = LSS.numOfFluids();
        continue;
      }
      tag[iNode] = 0;
      for(int i=0; i<dim; i++) {
        if(phi[iNode][i]>0.0) {
          if (programmedBurn && (programmedBurn->isUnburnedEOS(i+1,burnTag) ||
                                 programmedBurn->isBurnedEOS(i+1,burnTag)) ) {
            if (/*programmedBurn->isIgnited(burnTag) &&*/
                (programmedBurn->nodeInside(burnTag,iSub,iNode) || programmedBurn->isFinished(burnTag)))
              tag[iNode] = programmedBurn->getBurnedEOS(burnTag);
            else {
              tag[iNode] = programmedBurn->getUnburnedEOS(burnTag);
            }
            break;
          }
          else{
            tag[iNode] = i+1; break;
          }
        }
      }
    }

    for(int iNode=0; iNode<Phi.subSize(iSub); iNode++)
      if(tag[iNode] == 2) fprintf(stderr," Caught Id = 2.\n");
  }
}

//------------------------------------------------------------------------------

template<int dim> /*this dim is actually dimLS*/
void FluidSelector::checkLSConsistency(DistSVec<double,dim> &Phi)
{
  int numLocSub = Phi.numLocSub();
  int iSub;
#pragma omp parallel for
  for(iSub=0; iSub<numLocSub; ++iSub) {
    double (*phi)[dim] = Phi.subData(iSub);
    int *tag           = fluidId->subData(iSub);
    for(int i=0; i<Phi.subSize(iSub); i++) {
      if(tag[i]==0) {
        if(phi[i][dim-1]>0.0) {
          fprintf(stderr,"BUG: Inconsistency between fluidId (%d) and phi (%e). numPhases = %d.\n", tag[i], phi[i][dim-1], numPhases);
          exit(-1);}}
      else if(tag[i]==numPhases) {
        if(fabs(phi[i][dim-1])>1.0e-10) {
          fprintf(stderr,"BUG: Inconsistency between fluidId (%d) and phi (%e). numPhases = %d.\n", tag[i], phi[i][dim-1], numPhases);
          exit(-1);}}
      else if (tag[i] == numPhases-1) {
        if(fabs(phi[i][dim-1]<=0.0)) {
          fprintf(stderr,"BUG: Inconsistency between fluidId (%d) and phi (%e). numPhases = %d.\n", tag[i], phi[i][dim-1], numPhases);
          exit(-1);}}
    }
  }
}

//------------------------------------------------------------------------------

