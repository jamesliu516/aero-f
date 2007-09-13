#include <DistBcData.h>

#include <IoData.h>
#include <VarFcn.h>
#include <BcData.h>
#include <SubDomain.h>
#include <Domain.h>

#include <math.h>

//------------------------------------------------------------------------------

template<int dim>
DistBcData<dim>::DistBcData(IoData &ioData, VarFcn *varFcn, Domain *domain) :
  Xdot(domain->getNodeDistInfo()), Temp(domain->getNodeDistInfo()),vf(varFcn),
  Ufarin(domain->getNodeDistInfo()), Ufarout(domain->getNodeDistInfo()),
  Uface(domain->getFaceDistInfo()), Unode(domain->getNodeDistInfo()),
  Uinletnode(domain->getInletNodeDistInfo()), rotInfo(ioData.rotations.rotationMap.dataMap)
{

  this->numLocSub = domain->getNumLocSub();
  this->subDomain = domain->getSubDomain();
  this->com = domain->getCommunicator();

  subBcData = new BcData<dim>*[this->numLocSub];
#pragma omp parallel for
  for (int iSub=0; iSub<this->numLocSub; ++iSub)
    subBcData[iSub] = new BcData<dim>(this->Uface(iSub), this->Unode(iSub), this->Uinletnode(iSub), this->Ufarin(iSub), this->Ufarout(iSub));


  if(vf->getType()==VarFcn::GAS || vf->getType()==VarFcn::GASINGAS)
    boundaryFluid = GAS;
  else
    boundaryFluid = TAIT;

  angles[0] = ioData.bc.inlet.alpha;
  angles[1] = ioData.bc.inlet.beta;

  gravityOn = false;
  gravity = 0.0;
  depth = 0.0;
  ngravity[0] = 0.0;
  ngravity[1] = 0.0;
  ngravity[2] = 0.0;
  if(ioData.bc.hydro.type == BcsHydroData::GRAVITY){
    gravityOn = true;
    gravity   = ioData.bc.hydro.gravity;
    depth     = ioData.bc.hydro.depth;
    ngravity[0] = cos(ioData.bc.hydro.alpha)*cos(ioData.bc.hydro.beta);
    ngravity[1] = cos(ioData.bc.hydro.alpha)*sin(ioData.bc.hydro.beta);
    ngravity[2] = sin(ioData.bc.hydro.alpha);
  }


  Xdot      = 0.0;
  Temp      = ioData.bc.wall.temperature;
  this->Ufarin     = 0.0;
  this->Ufarout    = 0.0;
  this->Unode      = 0.0;
  this->Uface      = 0.0;
  this->Uinletnode = 0.0;
  tref = ioData.ref.rv.time;
  vref = ioData.ref.rv.velocity;

}

//------------------------------------------------------------------------------

template<int dim>
DistBcData<dim>::~DistBcData()
{
 
  if (subBcData) {
#pragma omp parallel for
    for (int iSub=0; iSub<this->numLocSub; ++iSub)
      if (subBcData[iSub]) 
	delete subBcData[iSub];
    delete [] subBcData;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void DistBcData<dim>::finalize(VarFcn *varFcn, DistSVec<double,3> &X)
{

  varFcn->conservativeToPrimitive(this->Uin, Vin);
  varFcn->conservativeToPrimitive(this->Uout, Vout);

  double Pressure = varFcn->getPressure(Vin);
  com->printf(2, "Inlet:");
  int k;
  for (k=0; k<dim; ++k)
    com->printf(2, " %g", Vin[k]);
  com->printf(2, " %g", Pressure);
  com->printf(2, "\nOutlet:");
  for (k=0; k<dim; ++k)
    com->printf(2, " %g", Vout[k]);
  com->printf(2,"\n");
  if(Ub[0]>0.0){
    com->printf(2, "Bubble or Shocktube(false):");
    for(k=0; k<dim; k++)
      com->printf(2, " %g", Ub[k]);
    com->printf(2,"\n");
  }

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->assignFreeStreamValues2(this->Ufarin(iSub),
                                 this->Ufarout(iSub), this->Uface(iSub),
                                 this->Uinletnode(iSub));
  update(X);
}

//------------------------------------------------------------------------------

template<int dim>
void DistBcData<dim>::update(DistSVec<double,3> &X)  {

#pragma omp parallel for
  for (int iSub=0; iSub<this->numLocSub; ++iSub) {
    double (*xdot)[3] = Xdot.subData(iSub);
    double *temp = Temp.subData(iSub);
    double (*unode)[dim] = this->Unode.subData(iSub);
    int *rotOwn = this->subDomain[iSub]->getRotOwn();
    NodeSet &nodes = this->subDomain[iSub]->getNodes();
    //int count = 0;
    if (rotOwn)  { 
      for (int i=0; i<this->Unode.subSize(iSub); ++i) {
        if (rotOwn[i]>=0) {    // node belongs to a (potential) "rotating" surface
          map<int,RotationData *>::iterator it = rotInfo.find(rotOwn[i]);
          if(it != rotInfo.end()) { // the rotation data have been defined
	    if(it->second->infRadius == RotationData::TRUE) {
              double vel = it->second->omega / vref;
	      unode[i][1] = vel*it->second->nx;
	      unode[i][2] = vel*it->second->ny;
	      unode[i][3] = vel*it->second->nz;
	    } 
            else {
	      double ox = tref*it->second->omega*it->second->nx;
	      double oy = tref*it->second->omega*it->second->ny;
	      double oz = tref*it->second->omega*it->second->nz;
	      double xd = nodes[i][0] - it->second->x0;
	      double yd = nodes[i][1] - it->second->y0;
	      double zd = nodes[i][2] - it->second->z0;
	      unode[i][1] = oy*zd-oz*yd + xdot[i][0];
	      unode[i][2] = oz*xd-ox*zd + xdot[i][1];
	      unode[i][3] = ox*yd-oy*xd + xdot[i][2];
	    }
	  } 
          else  { // no rotation data -> use velocity from mesh motion if any
            unode[i][1] = xdot[i][0];
            unode[i][2] = xdot[i][1];
            unode[i][3] = xdot[i][2];
	  }
        }
        unode[i][4] = temp[i];
      }
    }
    else { // node does not belong to a "rotating" surface
      for (int i=0; i<this->Unode.subSize(iSub); ++i) {
        unode[i][1] = xdot[i][0];
        unode[i][2] = xdot[i][1];
        unode[i][3] = xdot[i][2];
        unode[i][4] = temp[i];
      }
    }
    //if(count) { fprintf(stderr," In DistBcData<dim>::update(): subd %3d has %6d 'rotating' nodes\n",this->subDomain[iSub]->getGlobSubNum(),count); fflush(stderr); }
#if defined(STRONG_INLET_BC)
    this->subDomain[iSub]->setNodeBcValue(Vin, this->Unode(iSub));
#endif
    this->subDomain[iSub]->computeFaceBcValue(this->Unode(iSub), this->Uface(iSub));
  }

  if ( this->gravity > 0.0 )
    updateFarField(X);
}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::updateFarField(DistSVec<double,3> &X)
{
// it does not matter what kind of simulation is done.
// we only need to know which fluid is assumed to be at the
//    boundary.
  if (this->boundaryFluid == this->GAS)
    updateFarFieldGas(X);

  else if(this->boundaryFluid == this->TAIT)
    updateFarFieldLiquid(X);

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->assignFreeStreamValues2(this->Ufarin(iSub),
                                 this->Ufarout(iSub), this->Uface(iSub),
                                 this->Uinletnode(iSub));
}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::updateFarFieldGas(DistSVec<double,3> &X)
{

  // flow properties
  double gam = this->vf->getGamma();
  double Pstiff = this->vf->getPressureConstant();

#pragma parallel omp for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double ptempin, ptempout, un, velin2, velout2;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      ptempin  = this->Vin[4] + this->Vin[0] *this->gravity*un;
      ptempout = this->Vout[4]+ this->Vout[0]*this->gravity*un;
      velin2  = this->Vin[1]*this->Vin[1]+this->Vin[2]*this->Vin[2]+this->Vin[3]*this->Vin[3];
      velout2 = this->Vout[1]*this->Vout[1]+this->Vout[2]*this->Vout[2]+this->Vout[3]*this->Vout[3];

      uin[inode][4] = 0.5*this->Vin[0]*velin2 + (ptempin+gam*Pstiff)/(gam-1.0);
      uout[inode][4] = 0.5*this->Vout[0]*velout2 + (ptempout+gam*Pstiff)/(gam-1.0);

    }

  }

}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::updateFarFieldLiquid(DistSVec<double,3> &X)
{

  //flow properties
  double a = this->vf->getAlphaWater();
  double b = this->vf->getBetaWater();
  double P = this->vf->getPrefWater();
  double c = this->vf->getCv();
  double coeff = (b-1.0)*this->gravity/(a*b);

#pragma parallel omp for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double rtempin, rtempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      rtempin  = pow(this->Vin[0], b-1.0)+coeff*(un-this->depth);
      rtempout = pow(this->Vout[0],b-1.0)+coeff*(un-this->depth);
      rtempin  = pow(rtempin, 1.0/(b-1.0));
      rtempout = pow(rtempout,1.0/(b-1.0));

      uin[inode][0] = rtempin;
      uin[inode][1] = rtempin*this->Vin[1];
      uin[inode][2] = rtempin*this->Vin[2];
      uin[inode][3] = rtempin*this->Vin[3];
      uin[inode][4] = rtempin*this->Uin[4]/this->Vin[0];

      uout[inode][0] = rtempout;
      uout[inode][1] = rtempout*this->Vout[1];
      uout[inode][2] = rtempout*this->Vout[2];
      uout[inode][3] = rtempout*this->Vout[3];
      uout[inode][4] = rtempout*this->Uout[4]/this->Vout[0];

    }

  }
}

//------------------------------------------------------------------------------
template<int dim>
void DistBcDataEuler<dim>::setBoundaryConditionsGas(IoData &iod,
                                DistSVec<double,3> &X)
{
/* case 1: no gravity (typically aero simulation)
 *   pressure and density are constants and imposed as are
 *   velocity is constant but imposed as Mach*speed of sound
 *
 * case 2: gravity and depth (typically hydro simulation)
 *   density is constant and imposed as is
 *   pressure is a function of the depth, 
 *          ie P = P_surface + rho*g*(h-z)
 *   velocity is constant and imposed as Mach*speed of sound at depth "depth"
 *
 */

  // flow properties
  double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
  double Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant/iod.ref.rv.pressure;
  double rhoin = iod.bc.inlet.density;
  double rhoout = iod.bc.outlet.density;
  double pressurein  = iod.bc.inlet.pressure + rhoin*this->gravity*this->depth;
  double pressureout = iod.bc.outlet.pressure + rhoout*this->gravity*this->depth;

  double velin2 = 0.0;
  double velout2 = 0.0;
  if(iod.bc.inlet.mach >= 0.0 && iod.bc.outlet.mach >= 0.0){
    velin2 = gam * (pressurein + Pstiff)*
      iod.bc.inlet.mach*iod.bc.inlet.mach / rhoin;
    velout2 = gam * (pressureout + Pstiff)*
      iod.bc.outlet.mach*iod.bc.outlet.mach / rhoout;
  }else if(iod.bc.inlet.velocity >= 0.0 && iod.bc.outlet.velocity >= 0.0 &&
           iod.bc.inlet.velocity+iod.bc.outlet.velocity>0.0){
    velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    velout2 = iod.bc.outlet.velocity*iod.bc.outlet.velocity;
  }else{
    this->com->fprintf(stdout, " no proper velocity or mach number specified\n");
    exit(1);
  }
  double velin = sqrt(velin2);
  double velout= sqrt(velout2);
  
// computation of boundary values "on average", ie at node of coordinates (0,0,0)
  this->Uin[0] = rhoin;
  this->Uin[1] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * cos(iod.bc.inlet.beta);
  this->Uin[2] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * sin(iod.bc.inlet.beta);
  this->Uin[3] = this->Uin[0] * velin * sin(iod.bc.inlet.alpha);
  this->Uin[4] = (pressurein+gam*Pstiff)/(gam - 1.0) + 0.5 * this->Uin[0] * velin2;

  this->Uout[0] = rhoout;
  this->Uout[1] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * cos(iod.bc.outlet.beta);
  this->Uout[2] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * sin(iod.bc.outlet.beta);
  this->Uout[3] = this->Uout[0] * velout * sin(iod.bc.outlet.alpha);
  this->Uout[4] = (pressureout+gam*Pstiff)/(gam - 1.0) + 0.5 * this->Uout[0] * velout2;

// computation for each node according to its depth
// this will be passed in DistTimeState to initialize simulation 
#pragma parallel omp for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double ptempin, ptempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      ptempin  = pressurein  + rhoin *this->gravity*un;
      ptempout = pressureout + rhoout*this->gravity*un;

      uin[inode][0] = this->Uin[0];
      uin[inode][1] = this->Uin[1];
      uin[inode][2] = this->Uin[2];
      uin[inode][3] = this->Uin[3];
      uin[inode][4] = 0.5*rhoin*velin2 + (ptempin+gam*Pstiff)/(gam-1.0); 

      uout[inode][0] = this->Uout[0];
      uout[inode][1] = this->Uout[1];
      uout[inode][2] = this->Uout[2];
      uout[inode][3] = this->Uout[3];
      uout[inode][4] = 0.5*rhoout*velout2 + (ptempout+gam*Pstiff)/(gam-1.0);

    }
  }
  for (int idim=0; idim<dim; idim++)
    this->Ub[idim] = 0.0;

}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::setBoundaryConditionsLiquid(IoData &iod, VarFcn *vf,
                                DistSVec<double,3> &X)
{
  // flow properties
  double a = vf->getAlphaWater();
  double b = vf->getBetaWater();
  double c = vf->getCv();
  double P = vf->getPrefWater();
  double coeff = -(b-1.0)*this->gravity/(a*b);
  double rhoin = iod.bc.inlet.density;
  double rhoout = iod.bc.outlet.density;

  double velin2 = 0.0;
  double velout2 = 0.0;
  if(iod.bc.inlet.mach >= 0.0 && iod.bc.outlet.mach >= 0.0){
    velin2 = iod.bc.inlet.mach*iod.bc.inlet.mach * a*b*pow(rhoin, b-1.0);
    velout2 = iod.bc.outlet.mach*iod.bc.outlet.mach * a*b*pow(rhoout, b-1.0);
  }else if (iod.bc.inlet.velocity >= 0.0 && iod.bc.outlet.velocity >= 0.0 &&
            iod.bc.inlet.velocity+iod.bc.outlet.velocity>0.0){
    velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    velout2 = iod.bc.outlet.velocity*iod.bc.outlet.velocity;
  }else{
    this->com->fprintf(stdout, " no proper velocity or mach number specified\n");
    exit(1);
  }
  double velin  = sqrt(velin2);
  double velout = sqrt(velout2);

// computation of boundary values "on average", ie at node of coordinates (0,0,0)
  this->Uin[0] = rhoin;
  this->Uin[1] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * cos(iod.bc.inlet.beta);
  this->Uin[2] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * sin(iod.bc.inlet.beta);
  this->Uin[3] = this->Uin[0] * velin * sin(iod.bc.inlet.alpha);
  this->Uin[4] = this->Uin[0]*(c*iod.bc.inlet.temperature + 0.5 * velin2);

  this->Uout[0] = rhoout;
  this->Uout[1] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * cos(iod.bc.outlet.beta);
  this->Uout[2] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * sin(iod.bc.outlet.beta);
  this->Uout[3] = this->Uout[0] * velout * sin(iod.bc.outlet.alpha);
  this->Uout[4] = this->Uout[0]*(c*iod.bc.outlet.temperature + 0.5 * velout2);


// computation for each node according to its depth
// this will be passed in DistTimeState to initialize simulation 
#pragma parallel omp for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double rtempin, rtempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      rtempin  = pow(rhoin ,b-1.0)+coeff*(un-this->depth);
      rtempout = pow(rhoout,b-1.0)+coeff*(un-this->depth);
      rtempin  = pow(rtempin, 1.0/(b-1.0));
      rtempout = pow(rtempout,1.0/(b-1.0));

      uin[inode][0] = rtempin;
      uin[inode][1] = rtempin*this->Uin[1]/this->Uin[0];
      uin[inode][2] = rtempin*this->Uin[2]/this->Uin[0];
      uin[inode][3] = rtempin*this->Uin[3]/this->Uin[0];
      uin[inode][4] = rtempin*this->Uin[4]/this->Uin[0];

      uout[inode][0] = rtempout;
      uout[inode][1] = rtempout*this->Uout[1]/this->Uout[0];
      uout[inode][2] = rtempout*this->Uout[2]/this->Uout[0];
      uout[inode][3] = rtempout*this->Uout[3]/this->Uout[0];
      uout[inode][4] = rtempout*this->Uout[4]/this->Uout[0];

    }

  }
  for (int idim=0; idim<dim; idim++)
    this->Ub[idim] = 0.0;

}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::setBoundaryConditionsGasGas(IoData &iod,
                                DistSVec<double,3> &X)
{

// this->Uin set up with inlet values and properties of fluidModel1 = GasModel1
  double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
  double Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant/iod.ref.rv.pressure;
  double rhoin = iod.bc.inlet.density;
  double rhoout = iod.bc.outlet.density;
  double pressurein = iod.bc.inlet.pressure + rhoin*this->gravity*this->depth;
  double pressureout = iod.bc.outlet.pressure+ rhoout*this->gravity*this->depth;

  double velin2 = 0.0;
  double velout2 = 0.0;
  if(iod.bc.inlet.mach >= 0.0 && iod.bc.outlet.mach >= 0.0){
    velin2 = gam * (pressurein+Pstiff) * 
      iod.bc.inlet.mach*iod.bc.inlet.mach / rhoin;
    velout2 = gam * (pressureout+ Pstiff)*
      iod.bc.outlet.mach*iod.bc.outlet.mach / rhoout;
  }else if (iod.bc.inlet.velocity>= 0.0 && iod.bc.outlet.velocity >= 0.0  &&
            iod.bc.inlet.velocity+iod.bc.outlet.velocity>0.0){
    velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    velout2 = iod.bc.outlet.velocity*iod.bc.outlet.velocity;
  }else{
    this->com->fprintf(stdout, " no proper velocity or mach number specified\n");
    exit(1);
  }
  double velin  = sqrt(velin2);
  double velout = sqrt(velout2);

  this->Uin[0] = rhoin;
  this->Uin[1] = this->Uin[0]*velin*cos(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta);
  this->Uin[2] = this->Uin[0]*velin*cos(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta);
  this->Uin[3] = this->Uin[0]*velin*sin(iod.bc.inlet.alpha);
  this->Uin[4] = (pressurein+gam*Pstiff)/(gam-1.0) + 0.5 * this->Uin[0] * velin2;

  this->Uout[0] = rhoout;
  this->Uout[1] = this->Uout[0]*velout*cos(iod.bc.outlet.alpha)*cos(iod.bc.outlet.beta);
  this->Uout[2] = this->Uout[0]*velout*cos(iod.bc.outlet.alpha)*sin(iod.bc.outlet.beta);
  this->Uout[3] = this->Uout[0]*velout*sin(iod.bc.outlet.alpha);
  this->Uout[4] = (pressureout+gam*Pstiff)/(gam-1.0) + 0.5 * this->Uout[0] * velout2;


#pragma parallel omp for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double ptempin, ptempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      ptempin  = pressurein + rhoin *this->gravity*un;
      ptempout = pressureout+ rhoout*this->gravity*un;

      uin[inode][0] = this->Uin[0];
      uin[inode][1] = this->Uin[1];
      uin[inode][2] = this->Uin[2];
      uin[inode][3] = this->Uin[3];
      uin[inode][4] = 0.5*this->Uin[0]*velin2 + (ptempin+gam*Pstiff)/(gam-1.0); 

      uout[inode][0] = this->Uout[0];
      uout[inode][1] = this->Uout[1];
      uout[inode][2] = this->Uout[2];
      uout[inode][3] = this->Uout[3];
      uout[inode][4] = 0.5*this->Uout[0]*velout2 + (ptempout+gam*Pstiff)/(gam-1.0);

    }
  }



  if(iod.mf.problem == MultiFluidData::BUBBLE){
    // for bubble type of computation
// Ub set up with bubble values and properties of fluidModel2 = GasModel2
// outlet values are used for angles of attack
    gam = iod.eqs.fluidModel2.gasModel.specificHeatRatio;
    Pstiff = iod.eqs.fluidModel2.gasModel.pressureConstant/iod.ref.rv.pressure;
    if(iod.mf.icd.s1.mach >= 0.0)
      velout2 = gam * (iod.mf.icd.s1.p + Pstiff) *
        iod.mf.icd.s1.mach*iod.mf.icd.s1.mach / iod.mf.icd.s1.rho;
    else
      velout2 = iod.mf.icd.s1.vel*iod.mf.icd.s1.vel;
    velout = sqrt(velout2);
    this->Ub[0] = iod.mf.icd.s1.rho;
    this->Ub[1] = this->Ub[0] * velout * cos(iod.bc.outlet.alpha) * cos(iod.bc.outlet.beta);
    this->Ub[2] = this->Ub[0] * velout * cos(iod.bc.outlet.alpha) * sin(iod.bc.outlet.beta);
    this->Ub[3] = this->Ub[0] * velout * sin(iod.bc.outlet.alpha);
    this->Ub[4] = (iod.mf.icd.s1.p + gam*Pstiff)/(gam - 1.0)+ 0.5 * this->Ub[0] * velout2;
  }
  else if(iod.mf.problem == MultiFluidData::SHOCKTUBE){

     // for shock tube type of computation
// Ub set up with inlet values and properties of fluidModel2 = GasModel2
// fluidModel1 is on the right/outlet 
// fluidModel2 is on the left/inlet
    fprintf(stdout, "\n\nSHOCKTUBE PROBLEM\n\n");
    gam = iod.eqs.fluidModel2.gasModel.specificHeatRatio;
    Pstiff = iod.eqs.fluidModel2.gasModel.pressureConstant/iod.ref.rv.pressure;
    if(iod.bc.inlet.mach >= 0.0){
      velin2 = gam * (iod.bc.inlet.pressure+Pstiff) *
        iod.bc.inlet.mach*iod.bc.inlet.mach / iod.bc.inlet.density;
    }else{
      velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    }
    velin = sqrt(velin2);
    // for initialization
    this->Ub[0] = iod.bc.inlet.density;
    this->Ub[1] = this->Ub[0] * velin * cos(iod.bc.inlet.alpha) * cos(iod.bc.inlet.beta);
    this->Ub[2] = this->Ub[0] * velin * cos(iod.bc.inlet.alpha) * sin(iod.bc.inlet.beta);
    this->Ub[3] = this->Ub[0] * velin * sin(iod.bc.inlet.alpha);
    this->Ub[4] = 1.0/(gam - 1.0) * (iod.bc.inlet.pressure+gam*Pstiff) + 0.5 * this->Ub[0] * velin2;
  
    // for boundary conditions of shocktube problems
    this->Uin[0] = this->Ub[0];
    this->Uin[1] = this->Ub[1];
    this->Uin[2] = this->Ub[2];
    this->Uin[3] = this->Ub[3];
    this->Uin[4] = this->Ub[4];
  
#pragma parallel omp for
    for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
      double (*uin)[dim]  = this->Ufarin.subData(iSub);
  
      for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
        uin[inode][0] = this->Uin[0];
        uin[inode][1] = this->Uin[1];
        uin[inode][2] = this->Uin[2];
        uin[inode][3] = this->Uin[3];
        uin[inode][4] = this->Uin[4];
      }
    }
    // End shocktube setup
  }

}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::setBoundaryConditionsLiquidLiquid(IoData &iod, VarFcn *vf,
                                DistSVec<double,3> &X)
{

  // flow properties
  double a = vf->getAlphaWater();
  double b = vf->getBetaWater();
  double c = vf->getCv();
  double P = vf->getPrefWater();
  double coeff = (b-1.0)*this->gravity/(a*b);
  double rhoin = iod.bc.inlet.density;
  double rhoout = iod.bc.outlet.density;

  double velin2 = 0.0;
  double velout2 = 0.0;
  if(iod.bc.inlet.mach >= 0.0 && iod.bc.outlet.mach >= 0.0){
    velin2 = iod.bc.inlet.mach*iod.bc.inlet.mach * a*b*pow(rhoin, b-1.0);
    velout2 = iod.bc.outlet.mach*iod.bc.outlet.mach * a*b*pow(rhoout, b-1.0);
  }else if (iod.bc.inlet.velocity>= 0.0 && iod.bc.outlet.velocity >= 0.0 &&
            iod.bc.inlet.velocity+iod.bc.outlet.velocity>0.0){
    velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    velout2 = iod.bc.outlet.velocity*iod.bc.outlet.velocity;
  }else{
    this->com->fprintf(stdout, " no proper velocity or mach number specified\n");
    exit(1);
  }
  double velin  = sqrt(velin2);
  double velout = sqrt(velout2);

// Uin and Uout set up with iod.bc.let.values and properties of fluidModel1 = LiquidModel1
// computation of boundary values "on average", ie at node of coordinates (0,0,0)
  this->Uin[0] = rhoin;
  this->Uin[1] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * cos(iod.bc.inlet.beta);
  this->Uin[2] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * sin(iod.bc.inlet.beta);
  this->Uin[3] = this->Uin[0] * velin * sin(iod.bc.inlet.alpha);
  this->Uin[4] = this->Uin[0]*(c*iod.bc.inlet.temperature + 0.5 * velin2);

  this->Uout[0] = rhoout;
  this->Uout[1] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * cos(iod.bc.outlet.beta);
  this->Uout[2] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * sin(iod.bc.outlet.beta);
  this->Uout[3] = this->Uout[0] * velout * sin(iod.bc.outlet.alpha);
  this->Uout[4] = this->Uout[0]*(c*iod.bc.outlet.temperature + 0.5 * velout2);


// computation for each node according to its depth
// this will be passed in DistTimeState to initialize simulation
#pragma parallel omp for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double rtempin, rtempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      rtempin  = pow(rhoin ,b-1.0)-coeff*(un-this->depth);
      rtempout = pow(rhoout,b-1.0)-coeff*(un-this->depth);
      rtempin  = pow(rtempin, 1.0/(b-1.0));
      rtempout = pow(rtempout,1.0/(b-1.0));

      uin[inode][0] = rtempin;
      uin[inode][1] = rtempin*this->Uin[1]/this->Uin[0];
      uin[inode][2] = rtempin*this->Uin[2]/this->Uin[0];
      uin[inode][3] = rtempin*this->Uin[3]/this->Uin[0];
      uin[inode][4] = rtempin*this->Uin[4]/this->Uin[0];

      uout[inode][0] = rtempout;
      uout[inode][1] = rtempout*this->Uout[1]/this->Uout[0];
      uout[inode][2] = rtempout*this->Uout[2]/this->Uout[0];
      uout[inode][3] = rtempout*this->Uout[3]/this->Uout[0];
      uout[inode][4] = rtempout*this->Uout[4]/this->Uout[0];
    }
  }


  a = vf->getAlphaWaterbis();
  b = vf->getBetaWaterbis();
  c = vf->getCvbis();

  if(iod.mf.problem == MultiFluidData::BUBBLE){
  // Ub set up with iod.bc.outlet.values and properties of fluidModel2 = LiquidModel2
  //for bubble type of computation
    if(iod.mf.icd.s1.mach >= 0.0)
      velout2 = iod.mf.icd.s1.mach*iod.mf.icd.s1.mach * a*b*pow(iod.mf.icd.s1.rho, b-1.0);
    else
      velout2 = iod.mf.icd.s1.vel*iod.mf.icd.s1.vel;
    velout = sqrt(velout2);
    this->Ub[0] = iod.mf.icd.s1.rho;
    this->Ub[1] = this->Ub[0] * velout * cos(iod.bc.outlet.alpha) * cos(iod.bc.outlet.beta);
    this->Ub[2] = this->Ub[0] * velout * cos(iod.bc.outlet.alpha) * sin(iod.bc.outlet.beta);
    this->Ub[3] = this->Ub[0] * velout * sin(iod.bc.outlet.alpha);
    this->Ub[4] = this->Ub[0]*(c*iod.mf.icd.s1.t + 0.5 * velout2);
  }
  else if(iod.mf.problem == MultiFluidData::SHOCKTUBE){

/////////////////////////////////////////////////////////////////////////////////////
//                for shock tube type of computation
// Ub set up with inlet values and properties of fluidModel2 = LiquidModel2
// fluidModel1 is on the right/outlet 
// fluidModel2 is on the left/inlet
    fprintf(stdout, "\n\nSHOCKTUBE PROBLEM\n\n");
    if(iod.bc.inlet.mach >= 0.0)
      velin2 = iod.bc.inlet.mach*iod.bc.inlet.mach* a*b*pow(iod.bc.inlet.density, b-1.0);
    else
      velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    velin = sqrt(velin2);
    // for initialization
    this->Ub[0] = iod.bc.inlet.density;
    this->Ub[1] = this->Ub[0] * velin * cos(iod.bc.inlet.alpha) * cos(iod.bc.inlet.beta);
    this->Ub[2] = this->Ub[0] * velin * cos(iod.bc.inlet.alpha) * sin(iod.bc.inlet.beta);
    this->Ub[3] = this->Ub[0] * velin * sin(iod.bc.inlet.alpha);
    this->Ub[4] = this->Ub[0]*(c*iod.bc.inlet.temperature + 0.5*velin2);
  
    // for boundary conditions in shocktube problems
    this->Uin[0] = this->Ub[0];
    this->Uin[1] = this->Ub[1];
    this->Uin[2] = this->Ub[2];
    this->Uin[3] = this->Ub[3];
    this->Uin[4] = this->Ub[4];
  
#pragma parallel omp for
    for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
      double (*uin)[dim]  = this->Ufarin.subData(iSub);
  
      for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
        uin[inode][0] = this->Uin[0];
        uin[inode][1] = this->Uin[1];
        uin[inode][2] = this->Uin[2];
        uin[inode][3] = this->Uin[3];
        uin[inode][4] = this->Uin[4];
      }
    }
  // End of shocktube setup
  }
  
}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::setBoundaryConditionsGasLiquid(IoData &iod, VarFcn *vf,
				DistSVec<double,3> &X)
{

  // flow properties
  double a = vf->getAlphaWater();
  double b = vf->getBetaWater();
  double c = vf->getCv();
  double P = vf->getPrefWater();
  double coeff = (b-1.0)*this->gravity/(a*b);
  double rhoin = iod.bc.inlet.density;
  double rhoout = iod.bc.outlet.density;

  double velin2 = 0.0;
  double velout2 = 0.0;
  if(iod.bc.inlet.mach >= 0.0 && iod.bc.outlet.mach >= 0.0){
    velin2 = iod.bc.inlet.mach*iod.bc.inlet.mach * a*b*pow(rhoin, b-1.0);
    velout2 = iod.bc.outlet.mach*iod.bc.outlet.mach * a*b*pow(rhoout, b-1.0);
  }else if (iod.bc.inlet.velocity >= 0.0 && iod.bc.outlet.velocity >= 0.0 &&
            iod.bc.inlet.velocity+iod.bc.outlet.velocity>0.0){
    velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    velout2 = iod.bc.outlet.velocity*iod.bc.outlet.velocity;
  }else{
    this->com->fprintf(stdout, " no proper velocity or mach number specified\n");
    exit(1);
  }
  double velin  = sqrt(velin2);
  double velout = sqrt(velout2);

// Uin and Uout set up with iod.bc.let.values and properties of fluidModel1 = LiquidModel1
// computation of boundary values "on average", ie at node of coordinates (0,0,0)
  this->Uin[0] = rhoin;
  this->Uin[1] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * cos(iod.bc.inlet.beta);
  this->Uin[2] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * sin(iod.bc.inlet.beta);
  this->Uin[3] = this->Uin[0] * velin * sin(iod.bc.inlet.alpha);
  this->Uin[4] = this->Uin[0]*(c*iod.bc.inlet.temperature + 0.5 * velin2);

  this->Uout[0] = rhoout;
  this->Uout[1] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * cos(iod.bc.outlet.beta);
  this->Uout[2] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * sin(iod.bc.outlet.beta);
  this->Uout[3] = this->Uout[0] * velout * sin(iod.bc.outlet.alpha);
  this->Uout[4] = this->Uout[0]*(c*iod.bc.outlet.temperature + 0.5 * velout2);


// computation for each node according to its depth
// this will be passed in DistTimeState to initialize simulation
#pragma parallel omp for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double rtempin, rtempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      rtempin  = pow(rhoin ,b-1.0)+coeff*(un-this->depth);
      rtempout = pow(rhoout,b-1.0)+coeff*(un-this->depth);
      rtempin  = pow(rtempin, 1.0/(b-1.0));
      rtempout = pow(rtempout,1.0/(b-1.0));

      uin[inode][0] = rtempin;
      uin[inode][1] = rtempin*this->Uin[1]/this->Uin[0];
      uin[inode][2] = rtempin*this->Uin[2]/this->Uin[0];
      uin[inode][3] = rtempin*this->Uin[3]/this->Uin[0];
      uin[inode][4] = rtempin*this->Uin[4]/this->Uin[0];

      uout[inode][0] = rtempout;
      uout[inode][1] = rtempout*this->Uout[1]/this->Uout[0];
      uout[inode][2] = rtempout*this->Uout[2]/this->Uout[0];
      uout[inode][3] = rtempout*this->Uout[3]/this->Uout[0];
      uout[inode][4] = rtempout*this->Uout[4]/this->Uout[0];
    }
  }


  if(iod.mf.problem == MultiFluidData::BUBBLE){
//        for bubble type of computation
// Ub set up with iod.bc.outlet.values and properties of fluidModel2 = LiquidModel2
    double gam = iod.eqs.fluidModel2.gasModel.specificHeatRatio;
    double Pstiff = iod.eqs.fluidModel2.gasModel.pressureConstant/iod.ref.rv.pressure;
    if(iod.mf.icd.s1.mach >= 0.0)
      velout2 = gam * (iod.mf.icd.s1.p + Pstiff) *
        iod.mf.icd.s1.mach*iod.mf.icd.s1.mach / iod.mf.icd.s1.rho;
    else
      velout2 = iod.mf.icd.s1.vel*iod.mf.icd.s1.vel;
    velout = sqrt(velout2);
    
    this->Ub[0] = iod.mf.icd.s1.rho;
    this->Ub[1] = this->Ub[0] * velout * cos(iod.bc.outlet.alpha) * cos(iod.bc.outlet.beta);
    this->Ub[2] = this->Ub[0] * velout * cos(iod.bc.outlet.alpha) * sin(iod.bc.outlet.beta);
    this->Ub[3] = this->Ub[0] * velout * sin(iod.bc.outlet.alpha);
    this->Ub[4] = (iod.mf.icd.s1.p+gam*Pstiff)/(gam - 1.0) + 0.5 * this->Ub[0] * velout2;
		
  }
  else if(iod.mf.problem == MultiFluidData::SHOCKTUBE){

/////////////////////////////////////////////////////////////////////////////////////
//                for shock tube type of computation
// Ub set up with inlet values and properties of fluidModel2 = LiquidModel2
// fluidModel1 = Liquid is on the right/outlet 
// fluidModel2 = Gas    is on the left/inlet
    fprintf(stdout, "\n\nSHOCKTUBE PROBLEM\n\n");
  
    double gam = iod.eqs.fluidModel2.gasModel.specificHeatRatio;
    double Pstiff = iod.eqs.fluidModel2.gasModel.pressureConstant/iod.ref.rv.pressure;
    if(iod.bc.inlet.mach >= 0.0)
      velin2 = gam * (iod.bc.inlet.pressure+Pstiff) *
        iod.bc.inlet.mach*iod.bc.inlet.mach / iod.bc.inlet.density;
    else
      velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    velin = sqrt(velin2);
    
  // for initialization
    this->Ub[0] = iod.bc.inlet.density;
    this->Ub[1] = this->Ub[0] * velin * cos(iod.bc.inlet.alpha) * cos(iod.bc.inlet.beta);
    this->Ub[2] = this->Ub[0] * velin * cos(iod.bc.inlet.alpha) * sin(iod.bc.inlet.beta);
    this->Ub[3] = this->Ub[0] * velin * sin(iod.bc.inlet.alpha);
    this->Ub[4] = 1.0/(gam - 1.0) * (iod.bc.inlet.pressure+gam*Pstiff) + 0.5 * this->Ub[0] * velin2;
  
  // for boundary conditions of shocktube problems
    this->Uin[0] = this->Ub[0];
    this->Uin[1] = this->Ub[1];
    this->Uin[2] = this->Ub[2];
    this->Uin[3] = this->Ub[3];
    this->Uin[4] = this->Ub[4];
  
#pragma parallel omp for
    for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
      double (*uin)[dim]  = this->Ufarin.subData(iSub);
  
      for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
        uin[inode][0] = this->Uin[0];
        uin[inode][1] = this->Uin[1];
        uin[inode][2] = this->Uin[2];
        uin[inode][3] = this->Uin[3];
        uin[inode][4] = this->Uin[4];
      }
    }
    fprintf(stdout, "Shocktube:\n");
    fprintf(stdout, "     Inlet:  %e %e %e\n", iod.bc.inlet.density, velin, iod.bc.inlet.pressure);
    fprintf(stdout, "     Outlet: %e %e %e\n", this->Uout[0], this->Uout[1]/this->Uout[0], P+a*pow(this->Uout[0],b));
    // End shocktube setup
  }
}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::setBoundaryConditionsLiquidGas(IoData &iod, VarFcn *vf,
				DistSVec<double,3> &X)
{

// this->Uin set up with inlet values and properties of fluidModel1 = GasModel1
  double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
  double Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant/iod.ref.rv.pressure;
  double rhoin = iod.bc.inlet.density;
  double rhoout = iod.bc.outlet.density;
  double pressurein = iod.bc.inlet.pressure + rhoin*this->gravity*this->depth;
  double pressureout = iod.bc.outlet.pressure+ rhoout*this->gravity*this->depth;

  double velin2 = 0.0;
  double velout2 = 0.0;
  if(iod.bc.inlet.mach >= 0.0 && iod.bc.outlet.mach >= 0.0){
    velin2 = gam * (pressurein+Pstiff) * 
      iod.bc.inlet.mach*iod.bc.inlet.mach / rhoin;
    velout2 = gam * (pressureout+ Pstiff)*
      iod.bc.outlet.mach*iod.bc.outlet.mach / rhoout;
  }else if (iod.bc.inlet.velocity >= 0.0 && iod.bc.outlet.velocity >= 0.0 &&
            iod.bc.inlet.velocity+iod.bc.outlet.velocity>0.0){
    velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    velout2 = iod.bc.outlet.velocity*iod.bc.outlet.velocity;
  }else{
    this->com->fprintf(stdout, " no proper velocity or mach number specified\n");
    exit(1);
  }
  double velin  = sqrt(velin2);
  double velout = sqrt(velout2);

  this->Uin[0] = rhoin;
  this->Uin[1] = this->Uin[0]*velin*cos(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta);
  this->Uin[2] = this->Uin[0]*velin*cos(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta);
  this->Uin[3] = this->Uin[0]*velin*sin(iod.bc.inlet.alpha);
  this->Uin[4] = (pressurein+gam*Pstiff)/(gam-1.0) + 0.5 * this->Uin[0] * velin2;

  this->Uout[0] = rhoout;
  this->Uout[1] = this->Uout[0]*velout*cos(iod.bc.outlet.alpha)*cos(iod.bc.outlet.beta);
  this->Uout[2] = this->Uout[0]*velout*cos(iod.bc.outlet.alpha)*sin(iod.bc.outlet.beta);
  this->Uout[3] = this->Uout[0]*velout*sin(iod.bc.outlet.alpha);

#pragma parallel omp for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double ptempin, ptempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      ptempin  = pressurein + rhoin *this->gravity*un;
      ptempout = pressureout+ rhoout*this->gravity*un;

      uin[inode][0] = this->Uin[0];
      uin[inode][1] = this->Uin[1];
      uin[inode][2] = this->Uin[2];
      uin[inode][3] = this->Uin[3];
      uin[inode][4] = 0.5*this->Uin[0]*velin2 + (ptempin+gam*Pstiff)/(gam-1.0); 

      uout[inode][0] = this->Uout[0];
      uout[inode][1] = this->Uout[1];
      uout[inode][2] = this->Uout[2];
      uout[inode][3] = this->Uout[3];
      uout[inode][4] = 0.5*this->Uout[0]*velout2 + (ptempout+gam*Pstiff)/(gam-1.0);

    }
  }

  if(iod.mf.problem == MultiFluidData::BUBBLE){
//        for bubble type of computation
// Ub set up with iod.bc.outlet.values and properties of fluidModel2 = LiquidModel2
    double a = vf->getAlphaWater();
    double b = vf->getBetaWater();
    double c = vf->getCv();
  
    if(iod.mf.icd.s1.mach >= 0.0)
      velout2 = iod.mf.icd.s1.mach*iod.mf.icd.s1.mach * a*b*pow(iod.mf.icd.s1.rho, b-1.0);
    else
      velout2 = iod.mf.icd.s1.vel*iod.mf.icd.s1.vel;
    velout = sqrt(velout2);
    this->Ub[0] = iod.mf.icd.s1.rho;
    this->Ub[1] = this->Ub[0] * velout * cos(iod.bc.outlet.alpha) * cos(iod.bc.outlet.beta);
    this->Ub[2] = this->Ub[0] * velout * cos(iod.bc.outlet.alpha) * sin(iod.bc.outlet.beta);
    this->Ub[3] = this->Ub[0] * velout * sin(iod.bc.outlet.alpha);
    this->Ub[4] = this->Ub[0]*(c*iod.mf.icd.s1.t + 0.5 * velout2);
	
  }
  else if(iod.mf.problem == MultiFluidData::SHOCKTUBE){
    this->com->printf(2, "simulation for shocktube not relevant\n"); 
    this->com->printf(2, "inverse FluidModel and FluidModel2\n");
    exit(1);
  }

}

//------------------------------------------------------------------------------

template<int dim>
DistBcDataEuler<dim>::DistBcDataEuler(IoData &iod, VarFcn *vf, Domain *dom, DistSVec<double,3> &X) : 
  DistBcData<dim>(iod, vf, dom)
{
  if (iod.eqs.numPhase == 1){
    if (iod.eqs.fluidModel.fluid == FluidModelData::GAS)
      setBoundaryConditionsGas(iod, X);
    else if(iod.eqs.fluidModel.fluid == FluidModelData::LIQUID)
      setBoundaryConditionsLiquid(iod, vf, X);
    
  }else if (iod.eqs.numPhase == 2){
    if (iod.eqs.fluidModel.fluid == FluidModelData::GAS &&
        iod.eqs.fluidModel2.fluid == FluidModelData::GAS)
      setBoundaryConditionsGasGas(iod, X);
    
    else if (iod.eqs.fluidModel.fluid == FluidModelData::LIQUID &&
        iod.eqs.fluidModel2.fluid == FluidModelData::LIQUID)
      setBoundaryConditionsLiquidLiquid(iod, vf, X);
    
    else if (iod.eqs.fluidModel.fluid == FluidModelData::LIQUID &&
        iod.eqs.fluidModel2.fluid == FluidModelData::GAS)
      setBoundaryConditionsGasLiquid(iod, vf, X);
    
    else if (iod.eqs.fluidModel.fluid == FluidModelData::GAS &&
        iod.eqs.fluidModel2.fluid == FluidModelData::LIQUID)
      setBoundaryConditionsLiquidGas(iod, vf, X);
  }

  if (dim == 5)
    this->finalize(vf, X);

}

//------------------------------------------------------------------------------

template<int dim>
DistBcDataSA<dim>::DistBcDataSA(IoData &iod, VarFcn *vf, Domain *dom, DistSVec<double,3> &X) : 
  DistBcDataEuler<dim>(iod, vf, dom, X)
{

  tmp = 0;
  vec2Pat = 0;

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION) {
    tmp = new DistSVec<double,2>(dom->getNodeDistInfo());
    vec2Pat = new CommPattern<double>(dom->getSubTopo(), this->com, CommPattern<double>::CopyOnSend);
#pragma omp parallel for
    for (int iSub = 0; iSub<this->numLocSub; ++iSub)
      this->subDomain[iSub]->setComLenNodes(2, *vec2Pat);
    vec2Pat->finalize();
  }

  this->Uin[5] = this->Uin[0] * iod.bc.inlet.nutilde;
  this->Uout[5] = this->Uout[0] * iod.bc.outlet.nutilde;

#pragma parallel omp for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      uin[inode][5] = this->Uin[5];
      uout[inode][5] = this->Uout[5];
    }
  }

  this->finalize(vf, X);

}

//------------------------------------------------------------------------------

template<int dim>
DistBcDataSA<dim>::~DistBcDataSA()
{

  if (tmp) delete tmp;
  if (vec2Pat) delete vec2Pat;

}

//------------------------------------------------------------------------------
// unode[i][5] contains mutilde = rho*nutilde (=0.0 by default)

template<int dim>
void DistBcDataSA<dim>::computeNodeValue(DistSVec<double,3> &X)
{

  if (tmp) {
    int iSub;
#pragma omp parallel for
    for (iSub=0; iSub<this->numLocSub; ++iSub) {
      this->subDomain[iSub]->computeNodeBcValue(X(iSub), this->Uface(iSub), (*tmp)(iSub));
      this->subDomain[iSub]->sndData(*vec2Pat, tmp->subData(iSub));
    }

    vec2Pat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < this->numLocSub; ++iSub) {
      this->subDomain[iSub]->addRcvData(*vec2Pat, tmp->subData(iSub));
      double (*t)[2] = tmp->subData(iSub);
      double (*unode)[dim] = this->Unode.subData(iSub);
      for (int i=0; i<this->Unode.subSize(iSub); ++i) {
	if (t[i][0] != 0.0) {
	  double w = 1.0 / t[i][0];
	  unode[i][5] = w * t[i][1];
	}
      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
DistBcDataKE<dim>::DistBcDataKE(IoData &iod, VarFcn *vf, Domain *dom, DistSVec<double,3> &X) : 
  DistBcDataEuler<dim>(iod, vf, dom, X)
{

  tmp = new DistSVec<double,3>(dom->getNodeDistInfo());
  vec3Pat = dom->getVec3DPat();

  this->Uin[5] = this->Uin[0] * iod.bc.inlet.kenergy;
  this->Uin[6] = this->Uin[0] * iod.bc.inlet.eps;
  this->Uout[5] = this->Uout[0] * iod.bc.outlet.kenergy;
  this->Uout[6] = this->Uout[0] * iod.bc.outlet.eps;

#pragma parallel omp for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      uin[inode][5] = this->Uin[5];
      uin[inode][6] = this->Uin[6];
      uout[inode][5] = this->Uout[5];
      uout[inode][6] = this->Uout[6];
    }
  }

  this->finalize(vf, X);

}

//------------------------------------------------------------------------------

template<int dim>
DistBcDataKE<dim>::~DistBcDataKE()
{

  if (tmp) delete tmp;

}
//------------------------------------------------------------------------------
// unode[i][5] contains k and unode[i][6] contains eps

template<int dim>
void DistBcDataKE<dim>::computeNodeValue(DistSVec<double,3> &X)
{

  int iSub;

#pragma omp parallel for
  for (iSub=0; iSub<this->numLocSub; ++iSub) {
    this->subDomain[iSub]->computeNodeBcValue(X(iSub), this->Uface(iSub), (*tmp)(iSub));
    this->subDomain[iSub]->sndData(*vec3Pat, tmp->subData(iSub));
  }

  vec3Pat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub) {
    this->subDomain[iSub]->addRcvData(*vec3Pat, tmp->subData(iSub));
    double (*t)[3] = tmp->subData(iSub);
    double (*unode)[dim] = this->Unode.subData(iSub);
    for (int i=0; i<this->Unode.subSize(iSub); ++i) {
      if (t[i][0] != 0.0) {
	double w = 1.0 / t[i][0];
	unode[i][5] = w * t[i][1];
	unode[i][6] = w * t[i][2];
      }
    }
  }

}

//------------------------------------------------------------------------------
