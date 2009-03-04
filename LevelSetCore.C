#include <LevelSet.h>



//-------------------------------------------------------------------------
LevelSet::LevelSet(IoData &iod, Domain *dom):
  Phin(dom->getNodeDistInfo()), Phinm1(dom->getNodeDistInfo()), Phinm2(dom->getNodeDistInfo()),
  Psi(dom->getNodeDistInfo()), dt(dom->getNodeDistInfo()), PsiRes(dom->getNodeDistInfo()),
  w(dom->getNodeDistInfo()), domain(dom), Phi0(dom->getNodeDistInfo()),
  Tag(dom->getNodeDistInfo())
{

  com = domain->getCommunicator();
  numLocSub = domain->getNumLocSub();

  data = new TimeData(iod);

  bandlevel = iod.mf.bandlevel;
  localtime = iod.mf.localtime;
  subIt = iod.mf.subIt;
  cfl_psi = iod.mf.cfl;
  typeTracking = iod.mf.typeTracking;
  copy = iod.mf.copy;
  conv_eps = iod.mf.eps;
  diff = bool(iod.mf.outputdiff);
  if(typeTracking == MultiFluidData::GRADIENT){
    fprintf(stdout, "***Warning: if reinitialization in band --> problem!\n ***         You need to reinitialize in the whole domain with this method\n");
  }

  lsgrad  = new DistNodalGrad<1,double>(iod,dom,2);

}
//-------------------------------------------------------------------------

LevelSet::~LevelSet()
{
  if(data) delete data;
  if(lsgrad) delete lsgrad;
}

//---------------------------------------------------------------------------------------------------------

void LevelSet::setupPhiVolumesInitialConditions(IoData &iod, DistVec<double> &Phi){

  // loop on all Volumes to setup Phi_0 
  if(!iod.volumes.volumeMap.dataMap.empty()){
    map<int, VolumeData *>::iterator it;
    for (it=iod.volumes.volumeMap.dataMap.begin(); it!=iod.volumes.volumeMap.dataMap.end();it++)
      if(it->second->type==VolumeData::FLUID)
        //each volume (it->first) is setup using Input variables 'volumeInitialConditions'
        //                                 and equation of state 'fluidModel'
        domain->setupPhiVolumesInitialConditions(it->first, Phi);
  }

}

//---------------------------------------------------------------------------------------------------------

void LevelSet::setupPhiMultiFluidInitialConditions(IoData &iod, DistSVec<double,3> &X, DistVec<double> &Phi){

  // Note that Phi was already initialized either to 1.0 everywhere
  // (if there were no volumes specified in the input file)
  // or to +1.0 where a volid was specified equal to 0 
  //    and to -1.0 where a volid was specified different from 0.
  // Therefore the only initialization
  // left to do is for the spheres and other geometric shapes.
  // For that, Phi must be multiplied by the signed distance
  // to that geometric shape.

  if(iod.mf.initialConditions.nspheres>0)
    for(int i=0; i<iod.mf.initialConditions.nspheres; i++)
      domain->setupPhiMultiFluidInitialConditionsSphere(*(iod.mf.initialConditions.sphere[i]),X,Phi);

  if(iod.mf.initialConditions.nplanes>0)
    domain->setupPhiMultiFluidInitialConditionsPlane(iod.mf.initialConditions.p1,X,Phi);

}

//---------------------------------------------------------------------------------------------------------
void LevelSet::update(DistVec<double> &Phi)
{

  if (data->use_nm2 && data->exist_nm1) {
    Phinm2 = Phinm1;
    data->exist_nm2 = true;
  }
  if (data->use_nm1) {
    Phinm1 = Phin;
    data->exist_nm1 = true;
  }
  Phin = Phi;
}

//---------------------------------------------------------------------------------------------------------
void LevelSet::writeToDisk(char *name)
{

  DistSVec<double,1> WritePhi(domain->getNodeDistInfo(), reinterpret_cast<double (*)[1]>(Phin.data()));
  domain->writeVectorToFile(name, 0, 0.0, WritePhi);

  if (data->use_nm1){
    DistSVec<double,1> WritePhi1(domain->getNodeDistInfo(), reinterpret_cast<double (*)[1]>(Phinm1.data()));
    domain->writeVectorToFile(name, 1, 0.0, WritePhi1);
  }

  if (data->use_nm2){
    DistSVec<double,1> WritePhi2(domain->getNodeDistInfo(), reinterpret_cast<double (*)[1]>(Phinm2.data()));
    domain->writeVectorToFile(name, 2, 0.0, WritePhi2);
  }

}
//---------------------------------------------------------------------------------------------------------
bool LevelSet::checkConvergencePsi(int iteration, double &res0)
{

  double eps = 1.0e-6;
  double res = sqrt(PsiRes*PsiRes);

  if(iteration == 0)
    res0 = res;

  double target = eps*res0;
  if(res < target || iteration >= subIt || res < 1.e-12){
    com->fprintf(stdout, "*** ReinitLS: SubIt = %d (init = %e, res = %e, target = %e)\n", iteration, res0, res, target);
    return true;
  }

  return false;

}
