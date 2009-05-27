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
#include <VectorSet.h>

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

// Included (MB)
  if (iod.problem.alltype == ProblemData::_STEADY_SENSITIVITY_ANALYSIS_) {
    dV = new DistSVec<double,dim>(domain->getNodeDistInfo());
  }
  else {
    dV = 0;
  }

  tmp2 = 0;
  vec2Pat = 0;
  smag = 0;
  wale = 0;  
  vms = 0;
  dles = 0;
  dvms = 0;
  spaceOp = 0;
  mutOmu = 0;
  Cs = 0;
  CsDvms = 0;
  CsDles = 0;
  forceGen = 0;

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
      dles = new DistDynamicLESTerm<dim>(varFcn, iod, domain);
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

  int order;
  order = spaceOp->getSpaceOrder();

  if(order == 1)
  {
   nodalForceWeights[0] = 1.0; 
   nodalForceWeights[1] = 1.0;
  }
  else if(order == 2)
  {
   nodalForceWeights[0] = 22.0;
   nodalForceWeights[1] = 7.0;
  }
  else
  {
   com->fprintf(stderr, "Problem in Setting Order : Order is neither 1 or 2 \n");
   exit(0);
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
  if (dvms) delete dvms;
  if (spaceOp) delete spaceOp;
  if (mutOmu) delete mutOmu;
  if (Cs) delete Cs;
  if (CsDvms) delete CsDvms;
  if (CsDles) delete CsDles;

}

//------------------------------------------------------------------------------
// the nodal force F is *** NOT *** assembled

template<int dim>
void PostOperator<dim>::computeNodalForce(DistSVec<double,3> &X, DistSVec<double,dim> &U, 
					  DistVec<double> &Pin, DistSVec<double,3> &F,
					  DistVec<double> *Phi)
{

  varFcn->conservativeToPrimitive(U,*V,Phi);

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeNodalForce(postFcn, (*bcData)(iSub), (*geoState)(iSub),
				       X(iSub), (*V)(iSub), Pin(iSub), F(iSub));
  }

}

//------------------------------------------------------------------------------
// the nodal force F is *** NOT *** assembled

// Included (MB)
template<int dim>
void PostOperator<dim>::computeDerivativeOfNodalForce(DistSVec<double,3> &X, DistSVec<double,3> &dX,
                                                                                                    DistSVec<double,dim> &U, DistSVec<double,dim> &dU,
                                                                                                    DistVec<double> &Pin, double dS[3], DistSVec<double,3> &dF)
{

//Remark: Error mesage for pointers
  if (dV == 0) {
    fprintf(stderr, "*** Error: Varible dV does not exist!\n");
    exit(1);
  }

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    varFcn->conservativeToPrimitive(U(iSub), (*V)(iSub));
    varFcn->conservativeToPrimitiveDerivative(U(iSub), dU(iSub), (*V)(iSub), (*dV)(iSub));
    subDomain[iSub]->computeDerivativeOfNodalForce(postFcn, (*bcData)(iSub), (*geoState)(iSub),
				       X(iSub), dX(iSub), (*V)(iSub), (*dV)(iSub), Pin(iSub), dS, dF(iSub));
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

// Included (MB)
template<int dim>
void PostOperator<dim>::computeDerivativeOfNodalHeatPower(DistSVec<double,3>& X, DistSVec<double,3>& dX, DistSVec<double,dim>& U, DistSVec<double,dim>& dU, double dS[3], DistVec<double>& dP)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    varFcn->conservativeToPrimitive(U(iSub), (*V)(iSub));
    varFcn->conservativeToPrimitiveDerivative(U(iSub), dU(iSub), (*V)(iSub), (*dV)(iSub));
    subDomain[iSub]->computeDerivativeOfNodalHeatPower(postFcn, (*bcData)(iSub), (*geoState)(iSub),
					   X(iSub), dX(iSub), (*V)(iSub), (*dV)(iSub), dS, dP(iSub));
  }

}

//------------------------------------------------------------------------------
// computes the non-dimensional forces and moments

template<int dim>
void PostOperator<dim>::computeForceAndMoment(Vec3D &x0, DistSVec<double,3> &X, 
					      DistSVec<double,dim> &U, 
                                              DistVec<double> *Phi, Vec3D *Fi, 
					      Vec3D *Mi, Vec3D *Fv, Vec3D *Mv, int hydro, 
                                              VecSet< DistSVec<double,3> > *mX, Vec<double> *genCF)
{

// Phi must be a null pointer for single-phase flow
// Phi points to a DistVec<double> for multi-phase flow
  int iSurf;
  for(iSurf = 0; iSurf < numSurf; ++iSurf) {
    Fi[iSurf] = 0.0;
    Mi[iSurf] = 0.0;
    Fv[iSurf] = 0.0;
    Mv[iSurf] = 0.0;
  }

  varFcn->conservativeToPrimitive(U, *V, Phi);

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
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

    if (mX) {
      SubVecSet<DistSVec<double,3>, SVec<double,3> > subMX(mX, iSub);
      subDomain[iSub]->computeForceAndMoment(surfOutMap, postFcn, (*bcData)(iSub), (*geoState)(iSub), 
					     X(iSub), (*V)(iSub), x0, fi, mi, fv, mv, hydro, &subMX, genCF);
    }
    else 
      subDomain[iSub]->computeForceAndMoment(surfOutMap, postFcn, (*bcData)(iSub), (*geoState)(iSub),
                                             X(iSub), (*V)(iSub), x0, fi, mi, fv, mv, hydro, 0, 0);

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
  double F[3], M[3];
  if(forceGen != 0)
    forceGen->getForcesAndMoments(U, X, F, M);

  Fi[0][0] += F[0]; //assuming there is only one surface! (or the "embedded force" is associated to surf#1.)
  Fi[0][1] += F[1];
  Fi[0][2] += F[2];
  Mi[0][0] += M[0];
  Mi[0][1] += M[1];
  Mi[0][2] += M[2];

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

// Included (MB)
// computes the derivative of non-dimensional forces and moments

template<int dim>
void PostOperator<dim>::computeDerivativeOfForceAndMoment(Vec3D &x0, DistSVec<double,3> &X, DistSVec<double,3> &dX,
                                                                                                             DistSVec<double,dim> &U, DistSVec<double,dim> &dU, double dS[3],
                                                                                                             Vec3D *dFi, Vec3D *dMi, Vec3D *dFv, Vec3D *dMv, int hydro)
{

  int iSurf;
  for(iSurf = 0; iSurf < numSurf; ++iSurf) {
    dFi[iSurf] = 0.0;
    dMi[iSurf] = 0.0;
    dFv[iSurf] = 0.0;
    dMv[iSurf] = 0.0;
  }

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    varFcn->conservativeToPrimitive(U(iSub), (*V)(iSub));
    varFcn->conservativeToPrimitiveDerivative(U(iSub), dU(iSub), (*V)(iSub), (*dV)(iSub));
    Vec3D *dfi = new Vec3D[numSurf];
    Vec3D *dmi = new Vec3D[numSurf];
    Vec3D *dfv = new Vec3D[numSurf];
    Vec3D *dmv = new Vec3D[numSurf];
    for(iSurf = 0; iSurf < numSurf; ++iSurf) {
      dfi[iSurf] = 0.0;
      dmi[iSurf] = 0.0;
      dfv[iSurf] = 0.0;
      dmv[iSurf] = 0.0;
    }
    subDomain[iSub]->computeDerivativeOfForceAndMoment(surfOutMap, postFcn, (*bcData)(iSub), (*geoState)(iSub),
                                                       X(iSub), dX(iSub), (*V)(iSub), (*dV)(iSub), dS,
                                                       x0, dfi, dmi, dfv, dmv, hydro);

    for(iSurf = 0; iSurf < numSurf; ++iSurf) {
#pragma omp critical
      dFi[iSurf] += dfi[iSurf];
#pragma omp critical
      dMi[iSurf] += dmi[iSurf];
#pragma omp critical
      dFv[iSurf] += dfv[iSurf];
#pragma omp critical
      dMv[iSurf] += dmv[iSurf];
    }
    delete [] dfi;
    delete [] dmi;
    delete [] dfv;
    delete [] dmv;
  } 

  for(iSurf = 0; iSurf < numSurf; ++iSurf) {
#pragma omp critical
    double dCoef[12] = {dFi[iSurf][0], dFi[iSurf][1], dFi[iSurf][2],
                        dMi[iSurf][0], dMi[iSurf][1], dMi[iSurf][2],
		        dFv[iSurf][0], dFv[iSurf][1], dFv[iSurf][2],
		        dMv[iSurf][0], dMv[iSurf][1], dMv[iSurf][2]};
    com->globalSum(12, dCoef);

    dFi[iSurf][0] = dCoef[0]; 
    dFi[iSurf][1] = dCoef[1];
    dFi[iSurf][2] = dCoef[2];

    dMi[iSurf][0] = dCoef[3]; 
    dMi[iSurf][1] = dCoef[4];
    dMi[iSurf][2] = dCoef[5];

    dFv[iSurf][0] = dCoef[6]; 
    dFv[iSurf][1] = dCoef[7];
    dFv[iSurf][2] = dCoef[8];

    dMv[iSurf][0] = dCoef[9]; 
    dMv[iSurf][1] = dCoef[10];
    dMv[iSurf][2] = dCoef[11];

  }

  map<int, int>::iterator it;
  iSurf = 1;
  for (it = surfOutMap.begin(); it != surfOutMap.end(); it++)  {
    if (it->second > 0)  {
      dFi[0] += dFi[it->second];

      dMi[0] += dMi[it->second];

      dFv[0] += dFv[it->second];

      dMv[0] += dMv[it->second];
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

  if ((type == PostFcn::DELTA_PLUS) || (type == PostFcn::SKIN_FRICTION)) {
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
  } 

  else if (type == PostFcn::CSDLES) {
     if(!CsDles) CsDles = new DistVec<double>(domain->getNodeDistInfo());
     *CsDles = 0.0;
     varFcn->conservativeToPrimitive(U, *V);
     dles->computeCsValue(A, *bcData, X, *V, *CsDles);
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double* q = Q.subData(iSub);
      double* cs = (*CsDles).subData(iSub);
      for (int i=0; i<Q.subSize(iSub); ++i) {
        q[i]  = cs[i];
      }
    }
  } 

  else if (type == PostFcn::CSDVMS) {
    if(!CsDvms) CsDvms = new DistVec<double>(domain->getNodeDistInfo());
    *CsDvms = 0.0;
    spaceOp->computePostOpDVMS(X, A, U, CsDvms, timeState);
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double* q = Q.subData(iSub);
      double* cs = (*CsDvms).subData(iSub);
      for (int i=0; i<Q.subSize(iSub); ++i) {
        q[i]  = cs[i];
      }
    }
  }

  else if (type == PostFcn::MUT_OVER_MU) {
    if(!mutOmu) mutOmu = new DistVec<double>(domain->getNodeDistInfo());
    *mutOmu = 0.0;
    varFcn->conservativeToPrimitive(U, *V);
                                                                                                                          
    if(vms) {
      vms->computeMutOMu(A, X, *V, *mutOmu);
    }
    else if(smag) {
      domain->computeMutOMuSmag(smag, A, X, *V, *mutOmu);
    }
    else if(dles) {
      dles->computeMutOMu(A, *bcData, X, *V, *mutOmu);
    }
    else if(wale) {
       domain->computeMutOMuWale(wale, A, X, *V, *mutOmu);
    }
    else if(dvms) {
      if(!Cs) Cs = new DistVec<double>(domain->getNodeDistInfo());
      *Cs = 0.0;
      spaceOp->computePostOpDVMS(X, A, U, Cs, timeState);
      dvms->computeMutOMu(A, X, *V, *Cs, *mutOmu);
    }
    else {
       fprintf(stderr,"MuTOverMu option valid only for LES computations..  Aborting ....\n"); exit(1);	
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
void PostOperator<dim>::computeScalarQuantity(PostFcn::ScalarType type,
                                              DistSVec<double,3>& X,
                                              DistSVec<double,dim>& U,
                                              DistVec<double>& A,
                                              DistVec<double>& Q,
                                              DistTimeState<dim> *timeState,
                                              DistVec<double>& Phi)
{
  int iSub;

  if ((type == PostFcn::DELTA_PLUS) || (type == PostFcn::SKIN_FRICTION)) {
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
  }

  else if (type == PostFcn::CSDLES) {
     if(!CsDles) CsDles = new DistVec<double>(domain->getNodeDistInfo());
     *CsDles = 0.0;
     varFcn->conservativeToPrimitive(U, *V);
     dles->computeCsValue(A, *bcData, X, *V, *CsDles);
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double* q = Q.subData(iSub);
      double* cs = (*CsDles).subData(iSub);
      for (int i=0; i<Q.subSize(iSub); ++i) {
        q[i]  = cs[i];
      }
   }
  }

  else if (type == PostFcn::CSDVMS) {
    if(!CsDvms) CsDvms = new DistVec<double>(domain->getNodeDistInfo());
    *CsDvms = 0.0;
    spaceOp->computePostOpDVMS(X, A, U, CsDvms, timeState);
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double* q = Q.subData(iSub);
      double* cs = (*CsDvms).subData(iSub);
      for (int i=0; i<Q.subSize(iSub); ++i) {
        q[i]  = cs[i];
      }
    }
  }

  else if (type == PostFcn::MUT_OVER_MU) {
    if(!mutOmu) mutOmu = new DistVec<double>(domain->getNodeDistInfo());
    *mutOmu = 0.0;
    varFcn->conservativeToPrimitive(U, *V);


    if(vms) {
      vms->computeMutOMu(A, X, *V, *mutOmu);
    }
    else if(smag) {
      domain->computeMutOMuSmag(smag, A, X, *V, *mutOmu);
    }
    else if(dles) {
      dles->computeMutOMu(A, *bcData, X, *V, *mutOmu);
    }
    else if(wale) {
       domain->computeMutOMuWale(wale, A, X, *V, *mutOmu);
    }
    else if(dvms) {
      if(!Cs) Cs = new DistVec<double>(domain->getNodeDistInfo());
      *Cs = 0.0;
      spaceOp->computePostOpDVMS(X, A, U, Cs, timeState);
      dvms->computeMutOMu(A, X, *V, *Cs, *mutOmu);
    }
    else {
       fprintf(stderr,"MuTOverMu option valid only for LES computations..  Aborting ....\n"); exit(1);
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

  else if (type == PostFcn::PHILEVEL_STRUCTURE) {
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      subDomain[iSub]->computeNodeScalarQuantity(type, postFcn, (*V)(iSub), X(iSub), Q(iSub), Phi(iSub));
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

//---------------------------------------------------------------------------------
// Included (MB)
template<int dim>
void PostOperator<dim>::computeDerivativeOfScalarQuantity(PostFcn::ScalarDerivativeType type, double dS[3], DistSVec<double,3>& X, DistSVec<double,3>& dX, DistSVec<double,dim>& U, DistSVec<double,dim>& dU, DistVec<double>& dQ, DistTimeState<dim> *timeState)
{

  int iSub;

#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) {
    varFcn->conservativeToPrimitive(U(iSub), (*V)(iSub));
    varFcn->conservativeToPrimitiveDerivative(U(iSub), dU(iSub), (*V)(iSub), (*dV)(iSub));
    subDomain[iSub]->computeDerivativeOfNodeScalarQuantity(type, postFcn, dS, (*V)(iSub), (*dV)(iSub), X(iSub), dX(iSub), dQ(iSub));
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
                                                                                              
                                                                                              
  if ((type == PostFcn::DELTA_PLUS) || (type == PostFcn::SKIN_FRICTION)) {
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
                                                                                              
                                                                                              
  varFcn->conservativeToPrimitive(U, *V, &Phi);
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      //varFcn->conservativeToPrimitive(U(iSub), (*V)(iSub));
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
      double *phi = Phi.subData(iSub);
      for (int i=0; i<tmp3.subSize(iSub); ++i) {
        double v[dim];
        varFcn->conservativeToPrimitive(u[i], v, phi[i]);
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

// Included (MB)
template<int dim>
void PostOperator<dim>::computeDerivativeOfVectorQuantity(PostFcn::VectorDerivativeType type, DistSVec<double,3> &X, DistSVec<double,3> &dX, DistSVec<double,dim> &U, DistSVec<double,dim> &dU, DistSVec<double,3> &dQ)
{

  int iSub;

  if (type == PostFcn::DERIVATIVE_VELOCITY_VECTOR) {
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double (*u)[dim] = U.subData(iSub);
      double (*du)[dim] = dU.subData(iSub);
      double (*dq)[3] = dQ.subData(iSub);

      for (int i=0; i<dQ.subSize(iSub); ++i) {
	double v[dim];
        double dv[dim];
	varFcn->conservativeToPrimitive(u[i], v);
	varFcn->conservativeToPrimitiveDerivative(u[i], du[i], v, dv);
	Vec3D dVel = varFcn->getVelocity(dv);
	dq[i][0] = dVel[0];
	dq[i][1] = dVel[1];
	dq[i][2] = dVel[2];
      }
    }
  }
  else if (type == PostFcn::DERIVATIVE_DISPLACEMENT) {
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub)
      dQ(iSub) = dX(iSub);
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
                                            Vec3D &CMi, Vec3D &CFv, Vec3D &CMv, 
                                            VecSet< DistSVec<double,3> > *mX , DistVec<double> *genCF)  
{

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

// Included (MB)
template<int dim>
void PostOperator<dim>::checkVec(DistSVec<double,3>& V)  {

  for (int iSub=0; iSub<numLocSub; ++iSub) {
    subDomain[iSub]->checkVec(V(iSub));
  }

}

//------------------------------------------------------------------------------

