#include <PostOperator.h>

#include <IoData.h>
#include <VarFcn.h>
#include <DistBcData.h>
#include <DistGeoState.h>
#include <SubDomain.h>
#include <Domain.h>
#include <Vector3D.h>
#include <DistVector.h>
#include <SmagorinskyLESTerm.h>
#include <WaleLESTerm.h>
#include <DynamicLESTerm.h>
#include <DistDynamicLESTerm.h>
#include <DistDynamicVMSTerm.h>
#include <SpaceOperator.h>

//------------------------------------------------------------------------------

template<int dim>
PostOperator<dim>::PostOperator(IoData &iod, VarFcn *vf, DistBcData<dim> *bc, 
				DistGeoState *gs, Domain *dom, DistSVec<double,dim> *v) : 
  varFcn(vf), bcData(bc), geoState(gs), domain(dom)
{

  threshold = iod.schemes.ns.eps;
  pressInfty = iod.aero.pressure;
  refLengthSq = iod.ref.length * iod.ref.length;
  numLocSub = dom->getNumLocSub();
  subDomain = dom->getSubDomain();
  com = dom->getCommunicator();

  if (v) 
    V = v->alias();
  else 
    V = new DistSVec<double,dim>(dom->getNodeDistInfo());

  tmp2 = 0;
  vec2Pat = 0;
  smag = 0;
  wale = 0;  
  vms = 0;
  dles = 0;
  dlest = 0;
  dvms = 0;
  spaceOp = 0;

  spaceOp = new SpaceOperator<dim>(iod, vf, bc, gs, dom);

  if (iod.eqs.type == EquationsData::NAVIER_STOKES &&
    iod.eqs.tc.type == TurbulenceClosureData::LES) {
    if (iod.eqs.tc.les.type == LESModelData::VMS) {
      vms = new DistVMSLESTerm<dim>(varFcn, iod, domain);
    }
    else if (iod.eqs.tc.les.type == LESModelData::SMAGORINSKY) {
      smag = new SmagorinskyLESTerm(iod, varFcn);
    }
    else if (iod.eqs.tc.les.type == LESModelData::WALE) {
       wale = new WaleLESTerm(iod, varFcn);
    }     
    else if (iod.eqs.tc.les.type == LESModelData::DYNAMIC){
      dles = new DistDynamicLESTerm<dim>(iod, domain);
      dlest = new DynamicLESTerm(iod,varFcn);
    }
    else if (iod.eqs.tc.les.type == LESModelData::DYNAMICVMS){
      dvms = new DistDynamicVMSTerm<dim>(varFcn, iod, domain);
    }
  }

  if (iod.eqs.type == EquationsData::EULER)
    postFcn = new PostFcnEuler(iod, varFcn);
  else if (iod.eqs.type == EquationsData::NAVIER_STOKES) {
    if (iod.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
      if (iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS)
	postFcn = new PostFcnSA(iod, varFcn);
      else if (iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES)
	postFcn = new PostFcnDES(iod, varFcn);
      else if (iod.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE)
	postFcn = new PostFcnKE(iod, varFcn);
    }
    else
      postFcn = new PostFcnNS(iod, varFcn);
  }

  numSurf = 1;
  map<int, SurfaceData *> &sMap = iod.surfaces.surfaceMap.dataMap;
  map<int, SurfaceData *>::iterator it;
  for(it = sMap.begin(); it != sMap.end(); ++it) {
    if(it->second->computeForces == SurfaceData::TRUE
       || (it->second->computeForces == SurfaceData::UNSPECIFIED 
          && it->second->forceResults == SurfaceData::YES) ){
      if(it->second->forceResults == SurfaceData::YES)
        surfOutMap[it->first] = numSurf++;
      else {
        surfOutMap[it->first] = 0;
      }
    } 
    else if(it->second->computeForces == SurfaceData::FALSE)
      surfOutMap[it->first] = -1;  // We do not want the force computation 
    else
      surfOutMap[it->first] = -2; // We want the default behavior
  
  }

}

//------------------------------------------------------------------------------

template<int dim>
PostOperator<dim>::~PostOperator()
{

  if (V) delete V;
  if (tmp2) delete tmp2;
  if (vec2Pat) delete vec2Pat;
  if (postFcn) delete postFcn;
  if (vms) delete vms;
  if (smag) delete smag;
  if (wale) delete wale;
  if (dles) delete dles;
  if (dlest) delete dlest;
  if (dvms) delete dvms;
  if (spaceOp) delete spaceOp;

}

//------------------------------------------------------------------------------
// the nodal force F is *** NOT *** assembled

template<int dim>
void PostOperator<dim>::computeNodalForce(DistSVec<double,3> &X, DistSVec<double,dim> &U, 
					  DistVec<double> &Pin, DistSVec<double,3> &F)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    varFcn->conservativeToPrimitive(U(iSub), (*V)(iSub));
    subDomain[iSub]->computeNodalForce(postFcn, (*bcData)(iSub), (*geoState)(iSub),
				       X(iSub), (*V)(iSub), Pin(iSub), F(iSub));
  }

}

//------------------------------------------------------------------------------
// the nodal heat power is *** NOT *** assembled

template<int dim>
void PostOperator<dim>::computeNodalHeatPower(DistSVec<double,3>& X, DistSVec<double,dim>& U, 
					      DistVec<double>& P)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    varFcn->conservativeToPrimitive(U(iSub), (*V)(iSub));
    subDomain[iSub]->computeNodalHeatPower(postFcn, (*bcData)(iSub), (*geoState)(iSub),
					   X(iSub), (*V)(iSub), P(iSub));
  }

}

//------------------------------------------------------------------------------
// computes the non-dimensional forces and moments

template<int dim>
void PostOperator<dim>::computeForceAndMoment(Vec3D &x0, DistSVec<double,3> &X, 
					      DistSVec<double,dim> &U, Vec3D *Fi, 
					      Vec3D *Mi, Vec3D *Fv, Vec3D *Mv, int hydro)
{

  int iSurf;
  for(iSurf = 0; iSurf < numSurf; ++iSurf) {
    Fi[iSurf] = 0.0;
    Mi[iSurf] = 0.0;
    Fv[iSurf] = 0.0;
    Mv[iSurf] = 0.0;
  }

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    varFcn->conservativeToPrimitive(U(iSub), (*V)(iSub));
    Vec3D *fi = new Vec3D[numSurf];
    Vec3D *mi = new Vec3D[numSurf];
    Vec3D *fv = new Vec3D[numSurf];
    Vec3D *mv = new Vec3D[numSurf];
    for(iSurf = 0; iSurf < numSurf; ++iSurf) {
      fi[iSurf] = 0.0;
      mi[iSurf] = 0.0;
      fv[iSurf] = 0.0;
      mv[iSurf] = 0.0;
    }
    subDomain[iSub]->computeForceAndMoment(surfOutMap, postFcn, (*bcData)(iSub), (*geoState)(iSub), 
					   X(iSub), (*V)(iSub), x0, fi, mi, fv, mv, hydro);
    for(iSurf = 0; iSurf < numSurf; ++iSurf) {
#pragma omp critical
      Fi[iSurf] += fi[iSurf];
#pragma omp critical
      Mi[iSurf] += mi[iSurf];
#pragma omp critical
      Fv[iSurf] += fv[iSurf];
#pragma omp critical
      Mv[iSurf] += mv[iSurf];
    }
    delete [] fi;
    delete [] mi;
    delete [] fv;
    delete [] mv;
  } 

  for(iSurf = 0; iSurf < numSurf; ++iSurf) {
#pragma omp critical
    double coef[12] = {Fi[iSurf][0], Fi[iSurf][1], Fi[iSurf][2],
                       Mi[iSurf][0], Mi[iSurf][1], Mi[iSurf][2],
		       Fv[iSurf][0], Fv[iSurf][1], Fv[iSurf][2],
		       Mv[iSurf][0], Mv[iSurf][1], Mv[iSurf][2]};
    com->globalSum(12, coef);

    Fi[iSurf][0] = coef[0]; 
    Fi[iSurf][1] = coef[1];
    Fi[iSurf][2] = coef[2];

    Mi[iSurf][0] = coef[3]; 
    Mi[iSurf][1] = coef[4];
    Mi[iSurf][2] = coef[5];

    Fv[iSurf][0] = coef[6]; 
    Fv[iSurf][1] = coef[7];
    Fv[iSurf][2] = coef[8];

    Mv[iSurf][0] = coef[9]; 
    Mv[iSurf][1] = coef[10];
    Mv[iSurf][2] = coef[11];
 
  }

  map<int, int>::iterator it;
  iSurf = 1;
  for (it = surfOutMap.begin(); it != surfOutMap.end(); it++)  {
    if (it->second > 0)  {
      Fi[0] += Fi[it->second];

      Mi[0] += Mi[it->second];

      Fv[0] += Fv[it->second];

      Mv[0] += Mv[it->second];
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
double PostOperator<dim>::computeInterfaceWork(DistSVec<double,3>& X, 
					       DistSVec<double,dim>& U,
					       DistVec<double> &Pin)
{

  double E = 0.0;

#pragma omp parallel for reduction(+: E)
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    varFcn->conservativeToPrimitive(U(iSub), (*V)(iSub));
    E += subDomain[iSub]->computeInterfaceWork(postFcn, (*bcData)(iSub), (*geoState)(iSub), 
					       X(iSub), (*V)(iSub), Pin(iSub));
  }

  com->globalSum(1, &E);

  return E;

}

//------------------------------------------------------------------------------

template<int dim>
void PostOperator<dim>::computeScalarQuantity(PostFcn::ScalarType type, 
					      DistSVec<double,3>& X, 
					      DistSVec<double,dim>& U, 
                                              DistVec<double>& A,
					      DistVec<double>& Q,
                                              DistTimeState<dim> *timeState)
{

  int iSub;

  if (type == PostFcn::DELTA_PLUS) {
    if (!tmp2)
      tmp2 = new DistSVec<double,2>(domain->getNodeDistInfo());
    if (!vec2Pat) {
      vec2Pat = new CommPattern<double>(domain->getSubTopo(), com, 
					CommPattern<double>::CopyOnSend);
#pragma omp parallel for
      for (iSub = 0; iSub<numLocSub; ++iSub)
	subDomain[iSub]->setComLenNodes(2, *vec2Pat);
      vec2Pat->finalize();
    }
      
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      varFcn->conservativeToPrimitive(U(iSub), (*V)(iSub));
      subDomain[iSub]->computeFaceScalarQuantity(type, postFcn, (*bcData)(iSub), 
						 (*geoState)(iSub), X(iSub), 
						 (*V)(iSub), (*tmp2)(iSub));
      subDomain[iSub]->sndData(*vec2Pat, tmp2->subData(iSub));
    }

    vec2Pat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->addRcvData(*vec2Pat, tmp2->subData(iSub));

      double (*t)[2] = tmp2->subData(iSub);
      double* q = Q.subData(iSub);

      for (int i=0; i<Q.subSize(iSub); ++i) {
	if (t[i][0] != 0.0)
	  q[i] = t[i][1] / t[i][0];
	else
	  q[i] = 0.0;
      }
    }
  } 
  
  else if (type == PostFcn::VORTICITY) {
    DistSVec<double,6> R(domain->getNodeDistInfo());
    DistSVec<double,3> ddx(domain->getNodeDistInfo());
    DistSVec<double,3> ddy(domain->getNodeDistInfo());
    DistSVec<double,3> ddz(domain->getNodeDistInfo());
    DistSVec<double,3> tmp3(domain->getNodeDistInfo());
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      subDomain[iSub]->computeWeightsLeastSquaresEdgePart(X(iSub), R(iSub));
      subDomain[iSub]->sndData(*(domain->getWeightPat()), R.subData(iSub));
    }
    domain->getWeightPat()->exchange();
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      subDomain[iSub]->addRcvData(*(domain->getWeightPat()), R.subData(iSub));
      subDomain[iSub]->computeWeightsLeastSquaresNodePart(R(iSub));
      double (*u)[dim] = U.subData(iSub);
      double (*t3)[3] = tmp3.subData(iSub);
      double (*x)[3] = X.subData(iSub);
      for (int i=0; i<tmp3.subSize(iSub); ++i) {
	double v[dim];
	varFcn->conservativeToPrimitive(u[i], v);
	t3[i][0] = v[1];
	t3[i][1] = v[2];
	t3[i][2] = v[3];
      }
      subDomain[iSub]->computeGradientsLeastSquares(X(iSub), R(iSub), tmp3(iSub),
						    ddx(iSub), ddy(iSub), ddz(iSub));
    }
    domain->assemble(domain->getVec3DPat(), ddx);
    domain->assemble(domain->getVec3DPat(), ddy);
    domain->assemble(domain->getVec3DPat(), ddz);
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double (*dudx)[3] = ddx.subData(iSub);
      double (*dudy)[3] = ddy.subData(iSub);
      double (*dudz)[3] = ddz.subData(iSub);
      double* q = Q.subData(iSub);
      for (int i=0; i<Q.subSize(iSub); ++i) {
	double w0 = dudy[i][0] - dudx[i][1];
	double w1 = dudx[i][2] - dudz[i][0];
	double w2 = dudz[i][1] - dudy[i][2];
	q[i] = sqrt(w0*w0 + w1*w1 + w2*w2);
      }
    }
  }
 
  else if (type == PostFcn::PSENSOR) {
    DistSVec<double,6> R(domain->getNodeDistInfo());
    DistSVec<double,dim> ddx(domain->getNodeDistInfo());
    DistSVec<double,dim> ddy(domain->getNodeDistInfo());
    DistSVec<double,dim> ddz(domain->getNodeDistInfo());
    DistSVec<double,3> tmp3(domain->getNodeDistInfo());
    domain->computeWeightsLeastSquares(X, R);
    domain->computeGradientsLeastSquares(X, R, *V, ddx, ddy, ddz);
    domain->computePressureSensor(threshold, X, *V, ddx, ddy, ddz, tmp3, Q);
  } else if (type == PostFcn::CSDLES) {
    DistSVec<double,2> *CsDeltaSq;
    DistVec<double> *Cs;
    DistVec<double> *VolSum;
    CsDeltaSq = new DistSVec<double,2>(domain->getNodeDistInfo());
    Cs = new DistVec<double>(domain->getNodeDistInfo());
    VolSum = new DistVec<double>(domain->getNodeDistInfo());
    *CsDeltaSq = 0.0; *Cs = 0.0; *VolSum = 0.0;
    varFcn->conservativeToPrimitive(U, *V);
    dles->computeTestFilterValues(*CsDeltaSq, *VolSum, X, *V);
    domain->computeDynamicLESTerm(dlest, *CsDeltaSq, X, *Cs, *VolSum); //function has been overloaded
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double* q = Q.subData(iSub);
      double* cs = (*Cs).subData(iSub);
      for (int i=0; i<Q.subSize(iSub); ++i) {
        q[i]  = cs[i];
      }
    }	
    delete (CsDeltaSq); delete (VolSum);
  } 

  else if (type == PostFcn::CSDVMS) {
    DistVec<double> *Cs;
    Cs = new DistVec<double>(domain->getNodeDistInfo());
    *Cs = 0.0;
    spaceOp->computePostOpDVMS(X, A, U, Cs, timeState);
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double* q = Q.subData(iSub);
      double* cs = (*Cs).subData(iSub);
      for (int i=0; i<Q.subSize(iSub); ++i) {
        q[i]  = cs[i];
      }
    }
  }
                                                                                                                          
  else if (type == PostFcn::MUT_OVER_MU) {
    DistVec<double> *mutOmu;
    mutOmu = new DistVec<double>(domain->getNodeDistInfo());
    *mutOmu = 0.0;
    varFcn->conservativeToPrimitive(U, *V);
                                                                                                                          
    if(vms) {
      fprintf(stderr,"MuTOverMu not yet implemented for VMS-LES..  Aborting ....\n"); exit(1);
     // vms->obtainMutOverMu(X,V,mutOmu);
    }
    else if(smag) {
      domain->computeMutOMuSmag(smag, A, X, *V, *mutOmu);
    }
    else if(wale) {
       domain->computeMutOMuWale(wale, A, X, *V, *mutOmu);
    } 
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double* q = Q.subData(iSub);
      double* mtOm = (*mutOmu).subData(iSub);
      for (int i=0; i<Q.subSize(iSub); ++i) {
        q[i]  = mtOm[i];
      }
    }
  }

  else {
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      varFcn->conservativeToPrimitive(U(iSub), (*V)(iSub));
      subDomain[iSub]->computeNodeScalarQuantity(type, postFcn, (*V)(iSub), X(iSub), Q(iSub));
    }
  }

}

//------------------------------------------------------------------------------
template<int dim>
void PostOperator<dim>::computeCP(DistSVec<double,3>& X, DistSVec<double,dim>& U, Vec3D &cp)  {

  DistVec<double> Q(domain->getNodeDistInfo());
  DistVec<double> XP(domain->getNodeDistInfo());
  DistVec<double> YP(domain->getNodeDistInfo());
  DistVec<double> ZP(domain->getNodeDistInfo());
  Q = 0.0;
  XP = 0.0;
  YP = 0.0;
  ZP = 0.0;
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    varFcn->conservativeToPrimitive(U(iSub), (*V)(iSub));
    subDomain[iSub]->computeNodeScalarQuantity(PostFcn::DIFFPRESSURE, postFcn, (*V)(iSub), X(iSub), Q(iSub));
    subDomain[iSub]->computeXP(postFcn, (*V)(iSub), X(iSub), XP(iSub), 0);
    subDomain[iSub]->computeXP(postFcn, (*V)(iSub), X(iSub), YP(iSub), 1);
    subDomain[iSub]->computeXP(postFcn, (*V)(iSub), X(iSub), ZP(iSub), 2);
  }

  double xp = XP.sum();
  double yp = YP.sum();
  double zp = ZP.sum();
  double p = Q.sum();

  if (p == 0)  {
    cp[0] = 0;
    cp[1] = 0;
    cp[2] = 0;
  } 
  else {
  cp[0] = xp/p;
  cp[1] = yp/p;
  cp[2] = zp/p;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void PostOperator<dim>::computeScalarQuantity(PostFcn::ScalarType type,
                                              DistSVec<double,3>& X,
                                              DistSVec<double,dim>& U,
                                              DistVec<double>& Q,
                                              DistVec<double>& Phi)
{
                                                                                              
                                                                                              
  int iSub;
                                                                                              
                                                                                              
  if (type == PostFcn::DELTA_PLUS) {
    if (!tmp2)
      tmp2 = new DistSVec<double,2>(domain->getNodeDistInfo());
    if (!vec2Pat) {
      vec2Pat = new CommPattern<double>(domain->getSubTopo(), com,
                                        CommPattern<double>::CopyOnSend);
#pragma omp parallel for
      for (iSub = 0; iSub<numLocSub; ++iSub)
        subDomain[iSub]->setComLenNodes(2, *vec2Pat);
      vec2Pat->finalize();
    }
                                                                                              
                                                                                              
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      varFcn->conservativeToPrimitive(U(iSub), (*V)(iSub));
      subDomain[iSub]->computeFaceScalarQuantity(type, postFcn, (*bcData)(iSub),
                                                 (*geoState)(iSub), X(iSub),
                                                 (*V)(iSub), (*tmp2)(iSub));
      subDomain[iSub]->sndData(*vec2Pat, tmp2->subData(iSub));
    }
                                                                                              
                                                                                              
    vec2Pat->exchange();
                                                                                              
                                                                                              
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->addRcvData(*vec2Pat, tmp2->subData(iSub));
                                                                                              
                                                                                              
      double (*t)[2] = tmp2->subData(iSub);
      double* q = Q.subData(iSub);
                                                                                              
                                                                                              
      for (int i=0; i<Q.subSize(iSub); ++i) {
        if (t[i][0] != 0.0)
          q[i] = t[i][1] / t[i][0];
        else
          q[i] = 0.0;
      }
    }
  } else if (type == PostFcn::VORTICITY) {
    DistSVec<double,6> R(domain->getNodeDistInfo());
    DistSVec<double,3> ddx(domain->getNodeDistInfo());
    DistSVec<double,3> ddy(domain->getNodeDistInfo());
    DistSVec<double,3> ddz(domain->getNodeDistInfo());
    DistSVec<double,3> tmp3(domain->getNodeDistInfo());
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      subDomain[iSub]->computeWeightsLeastSquaresEdgePart(X(iSub), R(iSub));
      subDomain[iSub]->sndData(*(domain->getWeightPat()), R.subData(iSub));
    }
    domain->getWeightPat()->exchange();
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      subDomain[iSub]->addRcvData(*(domain->getWeightPat()), R.subData(iSub));
      subDomain[iSub]->computeWeightsLeastSquaresNodePart(R(iSub));
      double (*u)[dim] = U.subData(iSub);
      double (*t3)[3] = tmp3.subData(iSub);
      double (*x)[3] = X.subData(iSub);
      for (int i=0; i<tmp3.subSize(iSub); ++i) {
        double v[dim];
        varFcn->conservativeToPrimitive(u[i], v);
        t3[i][0] = v[1];
        t3[i][1] = v[2];
        t3[i][2] = v[3];
      }
      subDomain[iSub]->computeGradientsLeastSquares(X(iSub), R(iSub), tmp3(iSub),
                                                    ddx(iSub), ddy(iSub), ddz(iSub));
    }
    domain->assemble(domain->getVec3DPat(), ddx);
    domain->assemble(domain->getVec3DPat(), ddy);
    domain->assemble(domain->getVec3DPat(), ddz);
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double (*dudx)[3] = ddx.subData(iSub);
      double (*dudy)[3] = ddy.subData(iSub);
      double (*dudz)[3] = ddz.subData(iSub);
      double* q = Q.subData(iSub);
      for (int i=0; i<Q.subSize(iSub); ++i) {
        double w0 = dudy[i][0] - dudx[i][1];
        double w1 = dudx[i][2] - dudz[i][0];
        double w2 = dudz[i][1] - dudy[i][2];
        q[i] = sqrt(w0*w0 + w1*w1 + w2*w2);
      }
    }
  } else if (type == PostFcn::PSENSOR) {
    DistSVec<double,6> R(domain->getNodeDistInfo());
    DistSVec<double,dim> ddx(domain->getNodeDistInfo());
    DistSVec<double,dim> ddy(domain->getNodeDistInfo());
    DistSVec<double,dim> ddz(domain->getNodeDistInfo());
    DistSVec<double,3> tmp3(domain->getNodeDistInfo());
    domain->computeWeightsLeastSquares(X, R);
    domain->computeGradientsLeastSquares(X, R, *V, ddx, ddy, ddz);
    domain->computePressureSensor(threshold, X, *V, ddx, ddy, ddz, tmp3, Q);
  } else {
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      varFcn->conservativeToPrimitive(U(iSub), (*V)(iSub), &Phi(iSub));
      subDomain[iSub]->computeNodeScalarQuantity(type, postFcn, (*V)(iSub), X(iSub), Q(iSub), Phi(iSub));
    }
  }
                                                                                              
                                                                                              
}

//------------------------------------------------------------------------------

template<int dim>
void PostOperator<dim>::computeVectorQuantity(PostFcn::VectorType type, 
					      DistSVec<double,3> &X, 
					      DistSVec<double,dim> &U, 
					      DistSVec<double,3> &Q)
{

  int iSub;

  if (type == PostFcn::VELOCITY) {
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double (*u)[dim] = U.subData(iSub);
      double (*q)[3] = Q.subData(iSub);

      for (int i=0; i<Q.subSize(iSub); ++i) {
	double v[dim];
	varFcn->conservativeToPrimitive(u[i], v);
	Vec3D vel = varFcn->getVelocity(v);
	q[i][0] = vel[0];
	q[i][1] = vel[1];
	q[i][2] = vel[2];
      }
    }
  }
  else if (type == PostFcn::DISPLACEMENT || type == PostFcn::FLIGHTDISPLACEMENT || type == PostFcn::LOCALFLIGHTDISPLACEMENT) {
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub)
      subDomain[iSub]->computeDisplacement(X(iSub), Q(iSub));
  }
}

//------------------------------------------------------------------------------
                                                                                                                                                                                             
template<int dim>
void PostOperator<dim>::computeVectorQuantity(PostFcn::VectorType type,
                                              DistSVec<double,3> &X,
                                              DistSVec<double,dim> &U,
                                              DistSVec<double,3> &Q,
                                              DistVec<double> &Phi)
{
                                                                                                                                                                                             
  int iSub;
                                                                                                                                                                                             
  if (type == PostFcn::VELOCITY) {
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double (*u)[dim] = U.subData(iSub);
      double (*q)[3] = Q.subData(iSub);
      double (*phi) = Phi.subData(iSub);
                                                                                                                                                                                             
      for (int i=0; i<Q.subSize(iSub); ++i) {
        double v[dim];
        varFcn->conservativeToPrimitive(u[i], v, phi[i]);
        Vec3D vel = varFcn->getVelocity(v);
        q[i][0] = vel[0];
        q[i][1] = vel[1];
        q[i][2] = vel[2];
      }
    }
  }
  else if (type == PostFcn::DISPLACEMENT) {
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub)
      subDomain[iSub]->computeDisplacement(X(iSub), Q(iSub));
  }
                                                                                              
}


//------------------------------------------------------------------------------

// the nodal force DF is *** NOT *** assembled

template<int dim>
void PostOperator<dim>::computeForceDerivs(DistSVec<double,3> &X, DistSVec<double, dim> &U,
                                      DistSVec<double,dim> &deltaU, Vec<double> &modalF,
                                      VecSet< DistSVec<double,3> > &mX )  {

  varFcn->conservativeToPrimitive(U, *V);
  domain->computeForceDerivs(varFcn, X, *V, deltaU, modalF, mX);

}

//------------------------------------------------------------------------------

template<int dim>
void PostOperator<dim>::computeForceCoefficients(Vec3D &x0, DistSVec<double,3> &X,
                                            DistSVec<double,dim> &U, Vec3D &CFi,
                                            Vec3D &CMi, Vec3D &CFv, Vec3D &CMv)  {


  CFi = 0.0;
  CMi = 0.0;
  CFv = 0.0;
  CMv = 0.0;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    Vec3D locCFi, locCMi, locCFv, locCMv;
    varFcn->conservativeToPrimitive(U(iSub), (*V)(iSub));
    subDomain[iSub]->computeForceCoefficients(postFcn, x0, (*geoState)(iSub), (*bcData)(iSub),
                                     X(iSub), (*V)(iSub), pressInfty, locCFi, locCMi, locCFv, locCMv);

#pragma omp critical
    CFi += locCFi;
#pragma omp critical
    CMi += locCMi;
#pragma omp critical
    CFv += locCFv;
#pragma omp critical
    CMv += locCMv;

  }

  com->globalSum(3, CFi.v);
  com->globalSum(3, CMi.v);
  com->globalSum(3, CFv.v);
  com->globalSum(3, CMv.v);

  CFi *= 2.0 * refLengthSq;
  CMi *= 2.0 * refLengthSq;
  CFv *= 2.0 * refLengthSq;
  CMv *= 2.0 * refLengthSq;

}

//------------------------------------------------------------------------------

