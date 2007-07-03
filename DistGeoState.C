#include <DistGeoState.h>

#include <IoData.h>
#include <TimeData.h>
#include <Domain.h>
#include <GeoData.h>
#include <GeoState.h>
#include <Vector3D.h>
#include <DistVector.h>
#include <Communicator.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//------------------------------------------------------------------------------

DistGeoState::DistGeoState(IoData &ioData, Domain *dom) : data(ioData), domain(dom)
{

  numLocSub = domain->getNumLocSub();
  com = domain->getCommunicator();
  lscale = ioData.ref.rv.tlength;
  oolscale = 1.0 / lscale;

  if (data.use_n) {
    Xn = new DistSVec<double,3>(domain->getNodeDistInfo());
    ctrlVol_n = new DistVec<double>(domain->getNodeDistInfo());
  }
  else {
    Xn = 0;
    ctrlVol_n = 0;
  }

  if (data.use_nm1) {
    Xnm1 = new DistSVec<double,3>(domain->getNodeDistInfo());
    ctrlVol_nm1 = new DistVec<double>(domain->getNodeDistInfo());
  }
  else {
    Xnm1 = 0;
    ctrlVol_nm1 = 0;
  }

  if (data.use_nm2) {
    Xnm2 = new DistSVec<double,3>(domain->getNodeDistInfo());
    ctrlVol_nm2 = new DistVec<double>(domain->getNodeDistInfo());
  }
  else {
    Xnm2 = 0;
    ctrlVol_nm2 = 0;
  }

  Xdot = new DistSVec<double,3>(domain->getNodeDistInfo());

  d2wall = new DistVec<double>(domain->getNodeDistInfo());
  if (ioData.input.d2wall[0] != 0) {
    char* name = new char[strlen(ioData.input.prefix) + strlen(ioData.input.d2wall) + 1];
    sprintf(name, "%s%s", ioData.input.prefix, ioData.input.d2wall);
    DistSVec<double,1> d2w(d2wall->info(), reinterpret_cast<double (*)[1]>(d2wall->data()));
    domain->readVectorFromFile(name, 0, 0, d2w, &oolscale);
    delete [] name;
  }
  else
    *d2wall = 0.0;

  edgeNorm = new DistVec<Vec3D>(domain->getEdgeDistInfo());
  faceNorm = new DistVec<Vec3D>(domain->getFaceDistInfo());
  edgeNormVel = new DistVec<double>(domain->getEdgeDistInfo());
  faceNormVel = new DistVec<double>(domain->getFaceDistInfo());
  edgeNorm_nm1 = 0;
  faceNorm_nm1 = 0;
  edgeNormVel_nm1 = 0;
  faceNormVel_nm1 = 0;
  edgeNorm_nm2 = 0;
  faceNorm_nm2 = 0;
  edgeNormVel_nm2 = 0;
  faceNormVel_nm2 = 0;
  
  inletNodeNorm = new DistVec<Vec3D>(domain->getInletNodeDistInfo());
  numFaceNeighb = new DistVec<int>(domain->getInletNodeDistInfo());


  if (data.typeNormals == ImplicitData::SECOND_ORDER_GCL) {
    edgeNorm_nm1 = new DistVec<Vec3D>(domain->getEdgeDistInfo());
    faceNorm_nm1 = new DistVec<Vec3D>(domain->getFaceDistInfo());
    edgeNormVel_nm1 = new DistVec<double>(domain->getEdgeDistInfo());
    faceNormVel_nm1 = new DistVec<double>(domain->getFaceDistInfo());
  } 
  else if (data.typeNormals == ImplicitData::SECOND_ORDER_EZGCL) {
    edgeNormVel_nm1 = new DistVec<double>(domain->getEdgeDistInfo());
    faceNormVel_nm1 = new DistVec<double>(domain->getFaceDistInfo());
  } 
  else if (data.typeNormals == ImplicitData::THIRD_ORDER_EZGCL) {
    edgeNormVel_nm1 = new DistVec<double>(domain->getEdgeDistInfo());
    faceNormVel_nm1 = new DistVec<double>(domain->getFaceDistInfo());
    edgeNormVel_nm2 = new DistVec<double>(domain->getEdgeDistInfo());
    faceNormVel_nm2 = new DistVec<double>(domain->getFaceDistInfo());
  } 

  subGeoState = new GeoState*[numLocSub];

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    subGeoState[iSub] = 0;
    double* d2w = d2wall->subData(iSub);
    for (int i=0; i<d2wall->subSize(iSub); ++i)
      d2w[i] += ioData.bc.wall.delta;
  }

}

//------------------------------------------------------------------------------

DistGeoState::~DistGeoState()
{

  if (Xn) delete Xn;
  if (Xnm1) delete Xnm1;
  if (Xnm2) delete Xnm2;
  if (Xdot) delete Xdot;

  if (ctrlVol_n) delete ctrlVol_n;
  if (ctrlVol_nm1) delete ctrlVol_nm1;
  if (ctrlVol_nm2) delete ctrlVol_nm2;

  if (d2wall) delete d2wall;

  if (edgeNorm) delete edgeNorm;
  if (faceNorm) delete faceNorm;
  if (edgeNormVel) delete edgeNormVel;
  if (faceNormVel) delete faceNormVel;
  if (edgeNorm_nm1) delete edgeNorm_nm1;
  if (faceNorm_nm1) delete faceNorm_nm1;
  if (edgeNormVel_nm1) delete edgeNormVel_nm1;
  if (faceNormVel_nm1) delete faceNormVel_nm1;
  if (edgeNorm_nm2) delete edgeNorm_nm2;
  if (faceNorm_nm2) delete faceNorm_nm2;
  if (edgeNormVel_nm2) delete edgeNormVel_nm2;
  if (faceNormVel_nm2) delete faceNormVel_nm2;
  
  if (inletNodeNorm) delete inletNodeNorm;
  if (numFaceNeighb) delete numFaceNeighb;
 

  if (subGeoState) {
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub)
      if (subGeoState[iSub]) delete subGeoState[iSub];
    
    delete [] subGeoState;
  }

}

//------------------------------------------------------------------------------

void DistGeoState::setup(char *name, TimeData &timeData,
			DistSVec<double,3> *X, DistVec<double> *ctrlVol)
{
  setup1(name, X, ctrlVol);
  setup2(timeData);
}

//------------------------------------------------------------------------------

void DistGeoState::setup1(char *name, DistSVec<double,3> *X, DistVec<double> *ctrlVol)
{

  if (!data.use_n) {
    Xn = X->alias();
    ctrlVol_n = ctrlVol->alias();
  } 
  if (!data.use_nm1) {
    Xnm1 = Xn->alias();
    ctrlVol_nm1 = ctrlVol_n->alias();
  } 
  if (!data.use_nm2) {
    Xnm2 = Xnm1->alias();
    ctrlVol_nm2 = ctrlVol_nm1->alias();
  } 

  bool read_n = false;
  bool read_nm1 = false;
  bool read_nm2 = false;

  // This effectively `restarts' the geometry of the mesh, recovering the
  // last 3 positions.
  if (name[0] != 0) {
    read_n = domain->readVectorFromFile(name, 0, 0, *Xn, &oolscale);
    if (data.use_nm1)
      read_nm1 = domain->readVectorFromFile(name, 1, 0, *Xnm1, &oolscale);
    if (data.use_nm2)
      read_nm2 = domain->readVectorFromFile(name, 2, 0, *Xnm2, &oolscale);
  }

  if (!read_n)
    domain->getReferenceMeshPosition(*Xn);
  if (data.use_nm1 && !read_nm1)
    *Xnm1 = *Xn;
  if (data.use_nm2 && !read_nm2)
    *Xnm2 = *Xnm1;

  data.config = 0;
    
  domain->computeControlVolumes(lscale, *Xn, *ctrlVol_n);
  if (data.use_nm1)
    domain->computeControlVolumes(lscale, *Xnm1, *ctrlVol_nm1);
  if (data.use_nm2)
    domain->computeControlVolumes(lscale, *Xnm2, *ctrlVol_nm2);

  *X = *Xn;
  *ctrlVol = *ctrlVol_n;

  com->printf(2, "Control volume statistics: min=%.3e, max=%.3e, total=%.3e\n", 
	      ctrlVol_n->min(), ctrlVol_n->max(), ctrlVol_n->sum());
}

//-----------------------------------------------------------------------------

void DistGeoState::setup2(TimeData &timeData)
{

  double oodtnm1 = 1.0 / timeData.dt_nm1;
  double oodtnm2 = 1.0 / timeData.dt_nm2;
  if (data.typeVelocities == ImplicitData::ZERO)
    *Xdot = 0.0;
  else
    *Xdot = oodtnm1 * (*Xn - *Xnm1);

  domain->computeNormalsGCL1(*Xnm1, *Xn, *Xdot, *edgeNorm, *edgeNormVel, *faceNorm, *faceNormVel);
  domain->computeInletNormals(*inletNodeNorm, *faceNorm, *numFaceNeighb);

  if (data.typeNormals == ImplicitData::SECOND_ORDER_GCL)
    domain->computeNormalsGCL1(*Xnm1, *Xn, *Xdot, *edgeNorm_nm1, *edgeNormVel_nm1, 
			       *faceNorm_nm1, *faceNormVel_nm1);
  else if (data.typeNormals == ImplicitData::SECOND_ORDER_EZGCL)
    domain->computeNormalsEZGCL1(oodtnm1, *Xnm1, *Xn, *edgeNorm, *edgeNormVel_nm1, 
				 *faceNorm, *faceNormVel_nm1);
  else if (data.typeNormals == ImplicitData::THIRD_ORDER_EZGCL) {
    domain->computeNormalsEZGCL1(oodtnm2, *Xnm2, *Xnm1, *edgeNorm, *edgeNormVel_nm2, 
				 *faceNorm, *faceNormVel_nm2);
    domain->computeNormalsEZGCL1(oodtnm1, *Xnm1, *Xn, *edgeNorm, *edgeNormVel_nm1, 
				 *faceNorm, *faceNormVel_nm1);
  }

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    if (!subGeoState[iSub])
      subGeoState[iSub] = new GeoState(data, (*ctrlVol_n)(iSub), (*ctrlVol_nm1)(iSub),
				       (*ctrlVol_nm2)(iSub), (*d2wall)(iSub), 
				       (*edgeNorm)(iSub), (*faceNorm)(iSub), 
				       (*edgeNormVel)(iSub), (*faceNormVel)(iSub),
				       (*inletNodeNorm)(iSub), (*numFaceNeighb)(iSub));

}

//------------------------------------------------------------------------------

void DistGeoState::compute(TimeData &timeData, DistSVec<double,3> &Xsdot, 
			   DistSVec<double,3> &X, DistVec<double> &ctrlVol)
{

  data.config += 1;
    
  domain->computeControlVolumes(lscale, X, ctrlVol);

  domain->computeVelocities(data.typeVelocities, timeData, Xsdot, *Xnm1, *Xn, X, *Xdot);

  if (data.typeNormals == ImplicitData::FIRST_ORDER_GCL || 
      data.typeNormals == ImplicitData::SECOND_ORDER_GCL) {
    domain->computeNormalsGCL1(*Xn, X, *Xdot, *edgeNorm, *edgeNormVel, 
			       *faceNorm, *faceNormVel);
    if (data.typeNormals == ImplicitData::SECOND_ORDER_GCL)
      domain->computeNormalsGCL2(timeData, *edgeNorm, *edgeNorm_nm1, *edgeNormVel, 
				 *edgeNormVel_nm1, *faceNorm, *faceNorm_nm1,
				 *faceNormVel, *faceNormVel_nm1);
  }
  else if (data.typeNormals == ImplicitData::FIRST_ORDER_EZGCL || 
	   data.typeNormals == ImplicitData::SECOND_ORDER_EZGCL ||
	   data.typeNormals == ImplicitData::THIRD_ORDER_EZGCL) {
    domain->computeNormalsEZGCL1(1.0/timeData.dt_n, *Xn, X, *edgeNorm, *edgeNormVel, 
				 *faceNorm, *faceNormVel);
    if (data.typeNormals == ImplicitData::SECOND_ORDER_EZGCL)
      domain->computeNormalsEZGCL2(timeData, *edgeNormVel, *edgeNormVel_nm1, 
				   *faceNormVel, *faceNormVel_nm1);
    else if (data.typeNormals == ImplicitData::THIRD_ORDER_EZGCL)
      domain->computeNormalsEZGCL3(timeData, *edgeNormVel, *edgeNormVel_nm1, *edgeNormVel_nm2,
				   *faceNormVel, *faceNormVel_nm1, *faceNormVel_nm2);
  }
  else if (data.typeNormals == ImplicitData::CURRENT_CFG) {
    domain->computeNormalsGCL1(*Xn, *Xn, *Xdot, *edgeNorm, *edgeNormVel, 
			       *faceNorm, *faceNormVel);
  }
  else if (data.typeNormals == ImplicitData::LATEST_CFG) {
    domain->computeNormalsGCL1(X, X, *Xdot, *edgeNorm, *edgeNormVel, 
			       *faceNorm, *faceNormVel);
  }
  
  domain->computeInletNormals(*inletNodeNorm, *faceNorm, *numFaceNeighb);

}

//------------------------------------------------------------------------------

void DistGeoState::interpolate(double dt, double dtLeft, 
			       DistSVec<double,3> &Xs, DistSVec<double,3> &X)
{

  double alpha = dt / (dt + dtLeft) - 1.0;

  X = Xs + alpha * (Xs - *Xn);

}

//------------------------------------------------------------------------------

void DistGeoState::update(DistSVec<double,3> &X, DistVec<double> &ctrlVol)
{

  if (data.use_nm2) {
    *Xnm2 = *Xnm1;
    *ctrlVol_nm2 = *ctrlVol_nm1;
  }
  if (data.use_nm1) {
    *Xnm1 = *Xn;
    *ctrlVol_nm1 = *ctrlVol_n;
  }
  if (data.use_n) {
    *Xn = X;
    *ctrlVol_n = ctrlVol;
  }

}

//------------------------------------------------------------------------------

void DistGeoState::writeToDisk(char *name)
{

  if (data.use_n)
    domain->writeVectorToFile(name, 0, 0.0, *Xn, &lscale);
  if (data.use_nm1)
    domain->writeVectorToFile(name, 1, 0.0, *Xnm1, &lscale);
  if (data.use_nm2) 
    domain->writeVectorToFile(name, 2, 0.0, *Xnm2, &lscale);

}

//------------------------------------------------------------------------------
