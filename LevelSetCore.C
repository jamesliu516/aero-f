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
  if(typeTracking == MultiFluidData::GRADIENT){
		fprintf(stdout, "***Warning: if reinitialization in band --> problem!\n ***         You need to reinitialize in the whole domain with this method\n");
	}

  lsgrad  = new DistNodalGrad<1,double>(iod,dom,2);

}

//-------------------------------------------------------------------------

LevelSet::~LevelSet()
{
  if(data) delete data;
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
