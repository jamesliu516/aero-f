#include <TimeState.h>

#include <TimeData.h>
#include <GeoState.h>
#include <Vector.h>
#include <GenMatrix.h>
#include <BcDef.h>

//------------------------------------------------------------------------------

template<int dim>
TimeState<dim>::TimeState(TimeData &_data, Vec<double> &_dt, Vec<double> &_idti, 
			  Vec<double> &_idtv, SVec<double,dim> &_Un,
			  SVec<double,dim> &_Unm1, SVec<double,dim> &_Unm2, 
			  SVec<double,dim> &_Rn) : 
  data(_data), dt(_dt), idti(_idti), idtv(_idtv), Un(_Un), Unm1(_Unm1), Unm2(_Unm2), Rn(_Rn)
{

}

//------------------------------------------------------------------------------
// Add Time Derivative Term, d(AW)/dt to the flux F
// If running Non-Modal, adds invA*d(AW)/dt to the flux invA*F

template<int dim>
void TimeState<dim>::add_dAW_dt(bool *nodeFlag, GeoState &geoState, 
					Vec<double> &ctrlVol, SVec<double,dim> &Q, 
					SVec<double,dim> &R)
{

  Vec<double>& ctrlVol_n = geoState.getCtrlVol_n();
  Vec<double>& ctrlVol_nm1 = geoState.getCtrlVol_nm1();
  Vec<double>& ctrlVol_nm2 = geoState.getCtrlVol_nm2();

  double c_np1, c_n, c_nm1, c_nm2;
  for (int i=0; i<dt.size(); ++i) {

    double invDt = 1.0 / dt[i];
    if (data.use_modal == true)  {
      c_np1 = data.alpha_np1 * ctrlVol[i];
      c_n   = data.alpha_n * ctrlVol_n[i];
      c_nm1 = data.alpha_nm1 * ctrlVol_nm1[i];
      c_nm2 = data.alpha_nm2 * ctrlVol_nm2[i];
    }
    else  {
      double invCtrlVol = 1.0 / ctrlVol[i];
      c_np1 = data.alpha_np1;
      c_n   = data.alpha_n * ctrlVol_n[i] * invCtrlVol;
      c_nm1 = data.alpha_nm1 * ctrlVol_nm1[i] * invCtrlVol;
      c_nm2 = data.alpha_nm2 * ctrlVol_nm2[i] * invCtrlVol;
    }

    for (int k=0; k<dim; ++k) {
      double dAWdt = invDt * (c_np1*Q[i][k] + c_n*Un[i][k] +
                            c_nm1*Unm1[i][k] + c_nm2*Unm2[i][k]);
      if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON)
        R[i][k] = dAWdt + 0.5 * (R[i][k] + Rn[i][k]);
      else
        R[i][k] += dAWdt;
    }
  }
}

//-//------------------------------------------------------------------------------
// Add Time Derivative Term, d(AW)/dt to the flux F
// If running Non-Modal, adds invA*d(AW)/dt to the flux invA*F

template<int dim>
void TimeState<dim>::add_dAW_dtLS(bool *nodeFlag, GeoState &geoState, 
					Vec<double> &ctrlVol, Vec<double> &Q, 
					Vec<double> &R, Vec<double> &Q1, Vec<double> &Q2)
{

  Vec<double>& ctrlVol_n = geoState.getCtrlVol_n();
  Vec<double>& ctrlVol_nm1 = geoState.getCtrlVol_nm1();
  Vec<double>& ctrlVol_nm2 = geoState.getCtrlVol_nm2();

  double c_np1, c_n, c_nm1, c_nm2;
  for (int i=0; i<dt.size(); ++i) {

    double invDt = 1.0 / dt[i];
    if (data.use_modal == true)  {
      c_np1 = data.alpha_np1 * ctrlVol[i];
      c_n   = data.alpha_n * ctrlVol_n[i];
      c_nm1 = data.alpha_nm1 * ctrlVol_nm1[i];
      c_nm2 = data.alpha_nm2 * ctrlVol_nm2[i];
    }
    else  {
      double invCtrlVol = 1.0 / ctrlVol[i];
      c_np1 = data.alpha_np1;
      c_n   = data.alpha_n * ctrlVol_n[i] * invCtrlVol;
      c_nm1 = data.alpha_nm1 * ctrlVol_nm1[i] * invCtrlVol;
      c_nm2 = data.alpha_nm2 * ctrlVol_nm2[i] * invCtrlVol;
    }

                                                                                                                      
      double dAWdt = invDt * (c_np1*Q[i] + c_n*Q1[i] +
                              c_nm1*Q2[i]);
      if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON)
        R[i] = dAWdt + 0.5 * (R[i] + Rn[i][1]);
      else
        R[i] += dAWdt;

  }
}
//-----------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToJacobianNoPrec(bool *nodeFlag, Vec<double> &ctrlVol, GenMat<Scalar,neq> &A,
                                   SVec<double,dim> &U, VarFcn *vf, int* nodeType)
{

  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;
  if (data.use_modal == true && data.use_freq == false) A *= 0.5;

  if(!nodeType){
    for (int i=0; i<dt.size(); ++i) 
      addToJacobianNoPrecLocal(i, ctrlVol[i], U, A);

  }else{
    for (int i=0; i<dt.size(); ++i) 
      if(!(nodeType[i]==BC_INLET_MOVING || nodeType[i]==BC_OUTLET_MOVING ||
           nodeType[i]==BC_INLET_FIXED || nodeType[i]==BC_OUTLET_FIXED) )
        addToJacobianNoPrecLocal(i, ctrlVol[i], U, A);
  }
}
//------------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToJacobianNoPrecLocal(int i, double vol, 
					SVec<double,dim> &U, GenMat<Scalar,neq> &A)
{
  double c_np1;
  if (data.use_modal == true)
    c_np1 = data.alpha_np1 * vol / dt[i];
  else
    c_np1 = data.alpha_np1 / dt[i];

  Scalar *Aii = A.getElem_ii(i);
  for (int k=0; k<neq; ++k)
    Aii[k + k*neq] += c_np1;

}
//------------------------------------------------------------------------------
//  This Part Performs the Low Mach Steady State Preconditioning
//  Reference: Preconditioning Methods for Low-Speed Flows
//             By Turkel et. al. (ICASE Publication)

template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToJacobianGasPrec(bool *nodeFlag, Vec<double> &ctrlVol, GenMat<Scalar,neq> &A,
                                   SVec<double,dim> &U, VarFcn *vf, double gam, 
				   double pstiff, double beta, double k1, double cmach,
				   Vec<double> &irey, int* nodeType)
{
  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;
  if (data.use_modal == true && data.use_freq == false) A *= 0.5;
                                                                                                                           
  double c_np1;
  if(!nodeType){
    for (int i=0; i<dt.size(); ++i) 
      addToJacobianGasPrecLocal(i,ctrlVol[i],gam,pstiff,beta,k1,cmach,irey[i],U,A);

  }else{
    for (int i=0; i<dt.size(); ++i) 
      if(!(nodeType[i]==BC_INLET_MOVING || nodeType[i]==BC_OUTLET_MOVING ||
           nodeType[i]==BC_INLET_FIXED || nodeType[i]==BC_OUTLET_FIXED) )
        addToJacobianGasPrecLocal(i,ctrlVol[i],gam,pstiff,beta,k1,cmach,irey[i],U,A);

  }
}
//------------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToJacobianGasPrecLocal(int i, double vol, double gam, 
				double pstiff, double beta, double k1, double cmach,
				double irey, SVec<double,dim> &U, GenMat<Scalar,neq> &A)
{
  double c_np1;
  if (data.use_modal == true)
    c_np1 = data.alpha_np1 * vol / dt[i];
  else
    c_np1 = data.alpha_np1 / dt[i];

  Scalar *Aii = A.getElem_ii(i);

  if(neq<5){		//turbulence model equation in segregated solver
    for (int k=0; k<neq; ++k)
      Aii[k + k*neq] += c_np1;
  }else{	//Navier-Stokes (part of segregated turb model or alone) or fully coupled

    double ro = Un[i][0];
    double invRho = 1.0/ro;
    double u  = Un[i][1] * invRho;
    double v  = Un[i][2] * invRho;
    double w  = Un[i][3] * invRho;
    double u2 = u*u;
    double v2 = v*v;
    double w2 = w*w;
    double q2 = u2 + v2 + w2;
    double gam1 = gam - 1.0;
    double p  = gam1 * (Un[i][4] - 0.5 * ro * q2) - gam*pstiff;
    double c2 = gam*(p+pstiff)/ro;
    double locMach = sqrt(q2/c2); //local Preconditioning (ARL)
    beta = fmax(k1*locMach, beta);
    beta = fmin((1.0+sqrt(irey))*beta,cmach);

    double beta2 =   beta * beta;
    double qhat2 = (q2 * gam1)/2.0;
 
    double nu = qhat2/c2;
    double mu = (1.0/beta2) - 1.0;

    double Pinv[5][5] = { {nu*mu + 1.0,  -u*mu*gam1/c2,      -v*mu*gam1/c2,        -w*mu*gam1/c2,       mu*gam1/c2   },
                          {u*nu*mu,     1.0 - u2*mu*gam1/c2, -u*v*mu*gam1/c2,      -u*w*mu*gam1/c2,     u*mu*gam1/c2 },
                          {v*nu*mu,     -u*v*mu*gam1/c2 ,    1.0 - v2*mu*gam1/c2,  -v*w*mu*gam1/c2,     v*mu*gam1/c2 },
                          {w*nu*mu,     -u*w*mu*gam1/c2 ,    -v*w*mu*gam1/c2,      1.0 - w2*mu*gam1/c2, w*mu*gam1/c2 },
     	  	{0.5*mu*(1.0+nu)*q2,    -u*mu*(1+nu), -v*mu*(1+nu), -w*mu*(1+nu), (1.0/beta2)+mu*nu } };

    for (int l=0; l<5; ++l)
      for (int m=0; m<5; ++m)
        Aii[l*neq+m] += c_np1*Pinv[l][m];


    //turbulence preconditioning
    if(neq==6){
      double t1 = Un[i][5] * invRho;
      double mup = mu*t1*gam1/c2;
      double Pt[6] = {mu*nu*t1, -mup*u, -mup*v, -mup*w, mup, 1.0};
      for (int k=0; k<6; k++)
        Aii[neq*(neq-1)+k] += c_np1*Pt[k];

    }else if(neq==7){
      double t1 = Un[i][5] * invRho;
      double t2 = Un[i][6] * invRho;
      double mup1 = mu*t1*gam1/c2;
      double mup2 = mu*t2*gam1/c2;
      double Pt[2][7] = { {mu*nu*t1, -mup1*u, -mup1*v, -mup1*w, mup1, 1.0, 0.0},
                          {mu*nu*t2, -mup2*u, -mup2*v, -mup2*w, mup2, 0.0, 1.0} };
      for (int k=0; k<7; k++){
        Aii[neq*(neq-2)+k] += c_np1*Pt[0][k];
        Aii[neq*(neq-1)+k] += c_np1*Pt[1][k];
      }
    }
  }
}
//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToJacobianLiquidPrec(bool *nodeFlag, Vec<double> &ctrlVol, GenMat<Scalar,neq> &A,
                                   SVec<double,dim> &U, VarFcn *vf, double beta, double k1, double cmach, 
				   Vec<double> &irey, int* nodeType)
{
  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;
  if (data.use_modal == true && data.use_freq == false) A *= 0.5;
                                                                                                                           
  if(!nodeType){
    for (int i=0; i<dt.size(); ++i) 
      addToJacobianLiquidPrecLocal(i,ctrlVol[i],vf,beta,k1,cmach,irey[i],U,A);

  }else{
    for (int i=0; i<dt.size(); ++i) 
      if(!(nodeType[i]==BC_INLET_MOVING || nodeType[i]==BC_OUTLET_MOVING ||
           nodeType[i]==BC_INLET_FIXED || nodeType[i]==BC_OUTLET_FIXED) )
        addToJacobianLiquidPrecLocal(i,ctrlVol[i],vf,beta,k1,cmach,irey[i],U,A);
  }
}
//------------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToJacobianLiquidPrecLocal(int i, double vol, VarFcn *vf,
				double beta, double k1, double cmach, double irey,
				SVec<double,dim> &U, GenMat<Scalar,neq> &A)
{
// ARL : turbulence preconditioning never tested ...
  double c_np1;
  if (data.use_modal == true)
    c_np1 = data.alpha_np1 * vol / dt[i];
  else
    c_np1 = data.alpha_np1 / dt[i];

  Scalar *Aii = A.getElem_ii(i);
  if(neq<5){            //turbulence model equation in segregated solver
    for (int k=0; k<neq; ++k)
      Aii[k + k*neq] += c_np1;
  }else{        //Navier-Stokes (part of segregated turb model or alone) or fully coupled
    double V[dim];
    vf->conservativeToPrimitive(Un[i],V);
    double e = vf->computeRhoEnergy(V)/V[0];
    double locMach = vf->computeMachNumber(V); //local Preconditioning (ARL)
    beta = fmax(k1*locMach, beta);
    beta = fmin((1.0+sqrt(irey))*beta, cmach);
    double oobeta2 = 1.0/(beta*beta);
    double oobeta2m1 = oobeta2 - 1.0;

    double Pinv[dim];
    for (int j=0; j<dim; j++)
      Pinv[j] = oobeta2m1*V[j];
    Pinv[0] = oobeta2;
    Pinv[4] = oobeta2m1*e;
    /* The preconditioning matrix is:
     * Pinv[dim][dim] = { { oobeta2,           0.0, 0.0, 0.0, 0.0 , 0.0, 0.0 },
     *                    {(oobeta2-1.0)*V[1], 1.0, 0.0, 0.0, 0.0 , 0.0, 0.0 },
     *                    {(oobeta2-1.0)*V[2], 0.0, 1.0, 0.0, 0.0 , 0.0, 0.0 },
     *                    {(oobeta2-1.0)*V[3], 0.0, 0.0, 1.0, 0.0 , 0.0, 0.0 },
     *                    {(oobeta2-1.0)*e,    0.0, 0.0, 0.0, 1.0 , 0.0, 0.0 },
     *                    {(oobeta2-1.0)*V[5], 0.0, 0.0, 0.0, 0.0 , 1.0, 0.0 },
     *                    {(oobeta2-1.0)*V[6], 0.0, 0.0, 0.0, 0.0 , 0.0, 1.0 } };
     * Take the first 5-by-5 matrix to get the Euler preconditioner
     *      the first 6-by-6 matrix to get the "SA"  preconditioner
     *      the whole 7-by-7 matrix to get the "k-e" preconditioner
     */
 
    for (int j=1; j<dim; j++)
      Aii[j + j*neq] += c_np1;
    for (int j=0; j<dim; j++)
      Aii[j*neq] += c_np1*Pinv[j];
  }
}  
//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToH1(bool *nodeFlag, Vec<double> &ctrlVol, GenMat<Scalar,neq> &A)
{

  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;

  if (data.use_modal == true && data.use_freq == false) A *= 0.5;

  double c_np1;
  for (int i = 0; i < dt.size(); ++i) {

    if (nodeFlag && !nodeFlag[i]) continue;

    if (data.use_modal == true)
      if (data.use_freq == true)
        c_np1 = data.alpha_np1 * ctrlVol[i];
      else
        c_np1 = data.alpha_np1 * ctrlVol[i] / dt[i];
    else
      c_np1 = data.alpha_np1 / dt[i];

    Scalar *Aii = A.getElem_ii(i);

    for (int k=0; k<neq; ++k)
      Aii[k + k*neq] += c_np1;

  }

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToH1(bool *nodeFlag, Vec<double> &ctrlVol,
                GenMat<Scalar,neq> &A, Scalar shift)
{

//  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;

  if (data.use_modal == true && data.use_freq == false) A *= 0.5;

  Scalar c_np1;
  for (int i=0; i<dt.size(); ++i) {

    if (data.use_modal == true)  {
      if (data.use_freq == true)
        c_np1 = shift * ctrlVol[i];
      else
        c_np1 = data.alpha_np1 * ctrlVol[i] / dt[i];
    }
    else
      c_np1 = data.alpha_np1 / dt[i];

    Scalar *Aii = A.getElem_ii(i);

    for (int k=0; k<neq; ++k) {
      Aii[k + k*neq] += c_np1;
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar>
void TimeState<dim>::addToH2(bool *nodeFlag, VarFcn *varFcn, Vec<double> &ctrlVol,
			     SVec<double,dim> &V, GenMat<Scalar,dim> &A)
{

  double dfdUi[dim*dim], dfdVi[dim*dim];

  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;

  double coef = data.alpha_np1;
  if (data.use_modal == true && data.use_freq == false) {
    A *= 2.0;
    coef *= 3.0;
  }
  
  double c_np1;
  for (int i=0; i<dt.size(); ++i) {

    if (nodeFlag && !nodeFlag[i]) continue;

    if (data.use_freq == true)
      c_np1 = data.alpha_np1 * ctrlVol[i];
    else
      c_np1 = coef * ctrlVol[i] / dt[i];

    int k;
    for (k=0; k<dim*dim; ++k) dfdUi[k] = 0.0;
    for (k=0; k<dim; ++k) dfdUi[k + k*dim] = c_np1;

    varFcn->postMultiplyBydUdV(V[i], dfdUi, dfdVi);
  
    Scalar *Aii = A.getElem_ii(i);

    for (k=0; k<dim*dim; ++k) Aii[k] += dfdVi[k];

  }

}

//------------------------------------------------------------------------------
template<int dim>
template<class Scalar>
void TimeState<dim>::addToH2(bool *nodeFlag, VarFcn *varFcn, Vec<double> &ctrlVol,
                             SVec<double,dim> &V, GenMat<Scalar,dim> &A, Scalar coefVol, double coefA)
{

  Scalar dfdUi[dim*dim], dfdVi[dim*dim];

  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;

  //double coef = 1;
//data.alpha_np1;
  if (data.use_modal == true && data.use_freq == false) {
    A *= coefA;
  }

  Scalar c_np1;
  for (int i=0; i<dt.size(); ++i) {

    if (nodeFlag && !nodeFlag[i]) continue;

    if (data.use_freq == true)
      c_np1 = coefVol * ctrlVol[i];
    else
      c_np1 = coefVol * ctrlVol[i] / dt[i];

    int k;
    for (k=0; k<dim*dim; ++k) dfdUi[k] = 0.0;
    for (k=0; k<dim; ++k) dfdUi[k + k*dim] = c_np1;

    varFcn->postMultiplyBydUdV(V[i], dfdUi, dfdVi);

    Scalar *Aii = A.getElem_ii(i);

    for (k=0; k<dim*dim; ++k) Aii[k] += dfdVi[k];

  }

}



//------------------------------------------------------------------------------

template<int dim>
template<class Scalar>
void TimeState<dim>::addToH2(bool *nodeFlag, VarFcn *varFcn,
                Vec<double> &ctrlVol, SVec<double,dim> &V,
                GenMat<Scalar,dim> &A, Scalar shift)
{

  Scalar dfdUi[dim*dim], dfdVi[dim*dim];

  //if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;

  if (data.use_modal && data.use_freq == false) A *= 0.5;

  for (int i=0; i<dt.size(); ++i) {

    if (nodeFlag && nodeFlag[i] == 0) continue;
      Scalar c_np1;
      if (data.use_freq == true)
        c_np1 = shift*ctrlVol[i];
      else
        c_np1 = data.alpha_np1 * ctrlVol[i] / dt[i];

    int k;
    for (k=0; k<dim*dim; ++k) dfdUi[k] = 0.0;
    for (k=0; k<dim; ++k) dfdUi[k + k*dim] = c_np1;

    varFcn->postMultiplyBydUdV(V[i], dfdUi, dfdVi);

    Scalar *Aii = A.getElem_ii(i);

    for (k=0; k<dim*dim; ++k) Aii[k] += dfdVi[k];

  }

}

//------------------------------------------------------------------------------
                                                                                                         
template<int dim>
template<class Scalar>
void TimeState<dim>::addToH2LS(bool *nodeFlag, VarFcn *varFcn, Vec<double> &ctrlVol,
                             SVec<double,dim> &V, GenMat<Scalar,1> &A)
{
                                                                                                         
  double dfdUi[1], dfdVi[1];
                                                                                                         
  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;
                                                                                                         
  if (data.use_modal == true && data.use_freq == false) A *= 0.5;
                                                                                                         
  double c_np1;
  for (int i=0; i<dt.size(); ++i) {
                                                                                                         
    if (nodeFlag && !nodeFlag[i]) continue;
                                                                                                         
    if (data.use_freq == true)
      c_np1 = data.alpha_np1 * ctrlVol[i];
    else
      c_np1 = data.alpha_np1 * ctrlVol[i] / dt[i];
                                                                                                         
    int k;
    for (k=0; k<1; ++k) dfdUi[k] = 0.0;
    for (k=0; k<1; ++k) dfdUi[k + k] = c_np1;
                                                                                                         
//    varFcn->postMultiplyBydUdV(V[i], dfdUi, dfdVi);
                                                                                                         
    Scalar *Aii = A.getElem_ii(i);
                                                                                                         
    for (k=0; k<1; ++k) Aii[k] += dfdUi[k];
                                                                                                         
  }
                                                                                                         
}
//------------------------------------------------------------------------------
template<int dim>
template<class Scalar>
void TimeState<dim>::addToH2Minus(bool *nodeFlag, VarFcn *varFcn, Vec<double> &ctrlVol,
                                  SVec<double,dim> &V, GenMat<Scalar,dim> &A)
{

  double dfdUi[dim*dim], dfdVi[dim*dim];

  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;

  if (data.use_modal == true && data.use_freq == false) A *= -0.5;

  double c_np1;
  for (int i=0; i<dt.size(); ++i) {

    if (nodeFlag && nodeFlag[i] == 0) continue;

    if (data.use_freq)
      c_np1 = -2.0*data.alpha_np1 * ctrlVol[i];
    else
      c_np1 = data.alpha_np1 * ctrlVol[i] / dt[i];

    int k;
    for (k=0; k<dim*dim; ++k) dfdUi[k] = 0.0;
    for (k=0; k<dim; ++k) dfdUi[k + k*dim] = c_np1;

    varFcn->postMultiplyBydUdV(V[i], dfdUi, dfdVi);

    Scalar *Aii = A.getElem_ii(i);

    for (k=0; k<dim*dim; ++k) Aii[k] += dfdVi[k];

  }

}

//------------------------------------------------------------------------------
                                                                                                                          
template<int dim>
void TimeState<dim>::get_dW_dt(bool *nodeFlag, GeoState &geoState,
                               Vec<double> &ctrlVol, SVec<double,dim> &Q,
                               SVec<double,dim> &R)
{
                                                                                                                          
  Vec<double>& ctrlVol_n = geoState.getCtrlVol_n();
  Vec<double>& ctrlVol_nm1 = geoState.getCtrlVol_nm1();
  Vec<double>& ctrlVol_nm2 = geoState.getCtrlVol_nm2();
                                                                                                                          
  double c_np1, c_n, c_nm1, c_nm2;
                                                                                                                          
  for (int i=0; i<dt.size(); ++i) {
    if (!nodeFlag[i]) {
      double invDt = 1.0 / dt[i];
      c_np1 = data.alpha_np1*ctrlVol[i];
      c_n   = data.alpha_n * ctrlVol_n[i];
      c_nm1 = data.alpha_nm1 * ctrlVol_nm1[i];
      c_nm2 = data.alpha_nm2 * ctrlVol_nm2[i];
                                                                                                                          
      for (int k=0; k<dim; ++k) {
        double dWdt = invDt * (c_np1*Q[i][k] + c_n*Un[i][k] +
                               c_nm1*Unm1[i][k] + c_nm2*Unm2[i][k]);
        R[i][k] += dWdt;
      }
    }
  }
                                                                                                                          
}
                                                                                                                          
//------------------------------------------------------------------------------
                                                                                                                          
template<int dim>
void TimeState<dim>::get_dWBar_dt(bool *nodeFlag, GeoState &geoState,
                                  Vec<double> &ctrlVol, SVec<double,dim> &QBar,
                                  SVec<double,dim> &UnBar, SVec<double,dim> &Unm1Bar,
                                  SVec<double,dim> &Unm2Bar, SVec<double,dim> &R)
{
                                                                                                                          
  Vec<double>& ctrlVol_n = geoState.getCtrlVol_n();
  Vec<double>& ctrlVol_nm1 = geoState.getCtrlVol_nm1();
  Vec<double>& ctrlVol_nm2 = geoState.getCtrlVol_nm2();
                                                                                                                          
  double c_np1, c_n, c_nm1, c_nm2;
                                                                                                                          
  for (int i=0; i<dt.size(); ++i) {
    if (!nodeFlag[i]) {
      double invDt = 1.0 / dt[i];
      c_np1 = data.alpha_np1*ctrlVol[i];
      c_n   = data.alpha_n * ctrlVol_n[i];
      c_nm1 = data.alpha_nm1 * ctrlVol_nm1[i];
      c_nm2 = data.alpha_nm2 * ctrlVol_nm2[i];
                                                                                                                          
      for (int k=0; k<dim; ++k) {
        double dWdt = invDt * (c_np1*QBar[i][k] + c_n*UnBar[i][k] +
                               c_nm1*Unm1Bar[i][k] + c_nm2*Unm2Bar[i][k]);
        R[i][k] += dWdt;
      }
    }
  }
                                                                                                                          
}
                                                                                                                          
//------------------------------------------------------------------------------
