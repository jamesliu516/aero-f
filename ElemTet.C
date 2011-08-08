#include <Elem.h>
#include <Face.h>

#include <FemEquationTerm.h>
#include <MacroCell.h>
#include <VMSLESTerm.h>
#include <DynamicVMSTerm.h>
#include <SmagorinskyLESTerm.h>
#include <WaleLESTerm.h>
#include <DynamicLESTerm.h>
#include <GenMatrix.h>
#include <cmath>
#include <GeoState.h>
#include <BasicGeometry.h>

//------------------------------------------------------------------------------
//--------------functions in ElemTet class
//------------------------------------------------------------------------------

template<int dim>
void ElemTet::computeGalerkinTerm(FemEquationTerm *fet, SVec<double,3> &X, 
				  Vec<double> &d2wall, SVec<double,dim> &V, 
				  SVec<double,dim> &R, Vec<GhostPoint<dim>*> *ghostPoints,LevelSetStructure *LSS)
{
  // In the case of an embedded simulation, check if the tetrahedra is actually active
  bool isTetInactive=true,isAtTheInterface=true;
  if(ghostPoints) // Then LSS is also a non null pointer. It has already been checked in Domain.C
    {
      for(int i=0;i<4;++i) 
	{
	  isTetInactive    = (isTetInactive && !(LSS->isActive(0,nodeNum(i))));
	  isAtTheInterface = (isAtTheInterface && LSS->isActive(0,nodeNum(i)));
	}
      if(isTetInactive) return;
    }
  isAtTheInterface = !isAtTheInterface;

  // Common Part
  double dp1dxj[4][3];
  double vol = computeGradientP1Function(X, dp1dxj);

  double d2w[4] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)],
                   d2wall[nodeNum(2)], d2wall[nodeNum(3)]};
  double *v[4]  = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)], V[nodeNum(3)]};

  // Embedded Case. Replace states in the structure by Ghost States
  if(ghostPoints && isAtTheInterface) // Then LSS is also a non null pointer. It has already been checked in Domain.C
    {
      GhostPoint<dim> *gp;
      for(int i=0;i<4;++i)
	{
	  if(!(LSS->isActive(0,nodeNum(i))))
	    {
	      gp   = ghostPoints->operator[](nodeNum(i));
	      v[i] = gp->getPrimitiveState();
	    }
	}
      /*
      if(isAtTheInterface)
	{
	  cout<<"Node1: "<<nodeNum(0)<<" "<<LSS->isActive(0,nodeNum(0))<<" "<<v[0][0]<<" "<<v[0][1]<<" "<<v[0][2]<<" "<<v[0][3]<<endl;
	  cout<<"Node2: "<<nodeNum(1)<<" "<<LSS->isActive(0,nodeNum(1))<<" "<<v[1][0]<<" "<<v[1][1]<<" "<<v[1][2]<<" "<<v[1][3]<<endl;
	  cout<<"Node3: "<<nodeNum(2)<<" "<<LSS->isActive(0,nodeNum(2))<<" "<<v[2][0]<<" "<<v[2][1]<<" "<<v[2][2]<<" "<<v[2][3]<<endl;
	  cout<<"Node4: "<<nodeNum(3)<<" "<<LSS->isActive(0,nodeNum(3))<<" "<<v[3][0]<<" "<<v[3][1]<<" "<<v[3][2]<<" "<<v[3][3]<<endl;
	  std::cin.get();
	}
      */
    }

  double r[3][dim], s[dim], pr[12];
  bool porousTermExists =  fet->computeVolumeTerm(dp1dxj, d2w, v, reinterpret_cast<double *>(r),
                                                  s, pr, vol, X, nodeNum(), volume_id);

  if(ghostPoints && isAtTheInterface) // We don't want to update States associated to ghost points
    {
      for (int j=0; j<4; ++j) {
	int idx = nodeNum(j);
	if(LSS->isActive(0,idx))
	  {
	    for (int k=0; k<dim; ++k)
	      {
		R[idx][k] += vol * ( (r[0][k] * dp1dxj[j][0] + r[1][k] * dp1dxj[j][1] +
				      r[2][k] * dp1dxj[j][2]) - fourth * s[k] );
	      }
	  }
      }
      if (porousTermExists) {
	for (int j=0; j<4; ++j) {
	  int idx = nodeNum(j);
	  if(LSS->isActive(0,idx))
	    {
	      for (int k=1; k<4; ++k)
		{
		  R[idx][k] += pr[3*j+k-1];
		}
	    }
	}
      }
    }
  else // All the states are updated
    {
      for (int j=0; j<4; ++j) {
	int idx = nodeNum(j);
	for (int k=0; k<dim; ++k)
	  {
	    R[idx][k] += vol * ( (r[0][k] * dp1dxj[j][0] + r[1][k] * dp1dxj[j][1] +
				  r[2][k] * dp1dxj[j][2]) - fourth * s[k] );
	  }
      }
      if (porousTermExists) {
	for (int j=0; j<4; ++j) {
	  int idx = nodeNum(j);
	  for (int k=1; k<4; ++k)
	    {
	      R[idx][k] += pr[3*j+k-1];
	    }
	}
      }
    }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void ElemTet::computeDerivativeOfGalerkinTerm(FemEquationTerm *fet, SVec<double,3> &X, SVec<double,3> &dX,
			      Vec<double> &d2wall, SVec<double,dim> &V, SVec<double,dim> &dV, double dMach,
			      SVec<double,dim> &dR)
{

  double dp1dxj[4][3];
  double vol = computeGradientP1Function(X, dp1dxj);

  double ddp1dxj[4][3];
  double dvol = computeDerivativeOfGradientP1Function(X, dX, ddp1dxj);

  double d2w[4] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)],
		   d2wall[nodeNum(2)], d2wall[nodeNum(3)]};
  double *v[4] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)], V[nodeNum(3)]};
  double *dv[4] = {dV[nodeNum(0)], dV[nodeNum(1)], dV[nodeNum(2)], dV[nodeNum(3)]};

  double r[3][dim], s[dim], pr[12];
  bool porousTermExists =  fet->computeVolumeTerm(dp1dxj, d2w, v, reinterpret_cast<double *>(r),
                                                  s, pr, vol, X, nodeNum(), volume_id);

  double dr[3][dim], ds[dim], dpr[12];
  fet->computeDerivativeOfVolumeTerm(dp1dxj, ddp1dxj, d2w, v, dv, dMach, reinterpret_cast<double *>(dr), ds, dpr, dvol, X, nodeNum(), volume_id);

  for (int j=0; j<4; ++j) {
    int idx = nodeNum(j);
    for (int k=0; k<dim; ++k)
      dR[idx][k] += dvol * ( (r[0][k] * dp1dxj[j][0] + r[1][k] * dp1dxj[j][1] +
			    r[2][k] * dp1dxj[j][2]) - fourth * s[k] ) + vol * ( (dr[0][k] * dp1dxj[j][0] + r[0][k] * ddp1dxj[j][0] + dr[1][k] * dp1dxj[j][1] + r[1][k] * ddp1dxj[j][1] +
			    dr[2][k] * dp1dxj[j][2] + r[2][k] * ddp1dxj[j][2]) - fourth * ds[k] );
  }

  if (porousTermExists) {
    for (int j=0; j<4; ++j) {
      int idx = nodeNum(j);
      for (int k=1; k<4; ++k)
        dR[idx][k] += dpr[3*j+k-1];
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void ElemTet::computeP1Avg(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test, SVec<double,6> &Sij_Test, 
                           Vec<double> &modS_Test, SVec<double,8> &Eng_Test, SVec<double,3> &X, 
                           SVec<double,dim> &V, double gam, double R)

{

  double dp1dxj[4][3], dudxj[4][3], u[4][3];
  double NCG;
  double Int;                         // stores intermediate values
  double gam1 = gam - 1.0;
  double vol = computeGradientP1Function(X, dp1dxj);

  int i,j,k,l;                        // incrementers for the loops
  int i_mom = 0;
  int i_eng = 0;
  int i_sij = -1;


  // Multiplication factor for averaging //

  NCG = vol/4.0;  

  // adding up  the flow variables computed at the cg to each node //
  
  for (i=0; i<dim; ++i){
    Int = 0.0;
    for (j=0; j<4; ++j){
      Int += V[nodeNum(j)][i];
    }
    for (j=0; j<4; ++j){
      VCap[nodeNum(j)][i] += (NCG*Int);
    }
  }

  // incrementing the counter that stores the volume sum //

  for (i=0; i<4; ++i)
    Mom_Test[nodeNum(i)][i_mom] += vol;


  // adding rho_u computed at the cg to each node //

  for (i=1; i<4; ++i){
    Int = 0.0;
    for (j=0; j<4; ++j){
      Int += V[nodeNum(j)][0]*V[nodeNum(j)][i];
    }
    ++i_mom;
    for (j=0; j<4; ++j){
      Mom_Test[nodeNum(j)][i_mom] += (NCG*Int);
    }
  }

  // adding rho_u_u computed at the cg to each node //

  for (i=1; i<4; ++i){
    for (k=i; k<4; ++k){
      Int = 0.0;
      for (j=0; j<4; ++j){
        Int += V[nodeNum(j)][0]*V[nodeNum(j)][i]*V[nodeNum(j)][k];
      }
      ++i_mom;
      for (j=0; j<4; ++j){
        Mom_Test[nodeNum(j)][i_mom] += (NCG*Int);
      }
    }
  }

  // adding rho_s_p computed at the cg to each node //

  // step -1 : getting the velocities in to u matrix

    u[0][0] = V[nodeNum(0)][1];
    u[0][1] = V[nodeNum(0)][2];
    u[0][2] = V[nodeNum(0)][3];

    u[1][0] = V[nodeNum(1)][1];
    u[1][1] = V[nodeNum(1)][2];
    u[1][2] = V[nodeNum(1)][3];

    u[2][0] = V[nodeNum(2)][1];
    u[2][1] = V[nodeNum(2)][2];
    u[2][2] = V[nodeNum(2)][3];

    u[3][0] = V[nodeNum(3)][1];
    u[3][1] = V[nodeNum(3)][2];
    u[3][2] = V[nodeNum(3)][3];

  
  // step -2 : compute velocity gradients

    dudxj[0][0] = dp1dxj[0][0]*u[0][0] + dp1dxj[1][0]*u[1][0] +
	        dp1dxj[2][0]*u[2][0] + dp1dxj[3][0]*u[3][0];

    dudxj[0][1] = dp1dxj[0][1]*u[0][0] + dp1dxj[1][1]*u[1][0] +
          dp1dxj[2][1]*u[2][0] + dp1dxj[3][1]*u[3][0];

    dudxj[0][2] = dp1dxj[0][2]*u[0][0] + dp1dxj[1][2]*u[1][0] +
	    dp1dxj[2][2]*u[2][0] + dp1dxj[3][2]*u[3][0];

    dudxj[1][0] = dp1dxj[0][0]*u[0][1] + dp1dxj[1][0]*u[1][1] +
           dp1dxj[2][0]*u[2][1] + dp1dxj[3][0]*u[3][1];

    dudxj[1][1] = dp1dxj[0][1]*u[0][1] + dp1dxj[1][1]*u[1][1] +
	        dp1dxj[2][1]*u[2][1] + dp1dxj[3][1]*u[3][1];

    dudxj[1][2] = dp1dxj[0][2]*u[0][1] + dp1dxj[1][2]*u[1][1] +
	          dp1dxj[2][2]*u[2][1] + dp1dxj[3][2]*u[3][1];

    dudxj[2][0] = dp1dxj[0][0]*u[0][2] + dp1dxj[1][0]*u[1][2] +
	    dp1dxj[2][0]*u[2][2] + dp1dxj[3][0]*u[3][2];

    dudxj[2][1] = dp1dxj[0][1]*u[0][2] + dp1dxj[1][1]*u[1][2] +
	      dp1dxj[2][1]*u[2][2] + dp1dxj[3][1]*u[3][2];

    dudxj[2][2] = dp1dxj[0][2]*u[0][2] + dp1dxj[1][2]*u[1][2] +
	        dp1dxj[2][2]*u[2][2] + dp1dxj[3][2]*u[3][2];

  // step -3 : compute |S| i.e. sqrt2S2

    double S[3][3];

    S[0][0] = dudxj[0][0];
    S[1][1] = dudxj[1][1];
    S[2][2] = dudxj[2][2];

    S[0][1] = 0.5 * (dudxj[0][1] + dudxj[1][0]);
    S[0][2] = 0.5 * (dudxj[0][2] + dudxj[2][0]);
    S[1][2] = 0.5 * (dudxj[1][2] + dudxj[2][1]);

    S[1][0] = S[0][1];
    S[2][0] = S[0][2];
    S[2][1] = S[1][2];

    double S2 = (S[0][0]*S[0][0] + S[0][1]*S[0][1] + S[0][2]*S[0][2] + S[1][0]*S[1][0] +
                 S[1][1]*S[1][1] + S[1][2]*S[1][2] + S[2][0]*S[2][0] + S[2][1]*S[2][1] +
                 S[2][2]*S[2][2]);

    double sqrt2S2 = sqrt(2.0 * S2);

  // step-4 : compute Pij
    
    double Pij[6];
    Pij[0] =  (2./3.) * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
    Pij[1] =  (2./3.) * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
    Pij[2] =  (2./3.) * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
    Pij[3] =  (dudxj[1][0] + dudxj[0][1]);
    Pij[4] =  (dudxj[2][0] + dudxj[0][2]);
    Pij[5] =  (dudxj[2][1] + dudxj[1][2]);

  // step-5: last step of assembling

  for (i=0; i<6; ++i){
    Int = 0.0;
    for (j=0; j<4; ++j){
      Int += V[nodeNum(j)][0]*sqrt2S2*Pij[i];
    }
    ++i_mom;
    for (j=0; j<4; ++j){
      Mom_Test[nodeNum(j)][i_mom] += (NCG*Int);
    }
  }  

  // computing filtered Sij  //

  for (i=0; i<6; ++i) {
    ++i_sij;
    for (j=0; j<4; ++j)
      Sij_Test[nodeNum(j)][i_sij] += vol*Pij[i];
  }
 
  // computing filtered |Sij| //

  for (j=0; j<4; ++j)
    modS_Test[nodeNum(j)] += vol*sqrt2S2;


  // adding rho_e computed at the cg to each node //


  double squ[4];

  for(i=0; i<4; ++i)
    squ[i] = u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2]; 
  
  Int = 0.0;
  for (j=0; j<4; ++j){
    Int += V[nodeNum(j)][0]*((1.0/gam1)*(V[nodeNum(j)][4]/V[nodeNum(j)][0])+0.5*squ[j]);
  }
  for (j=0; j<4; ++j){
    Eng_Test[nodeNum(j)][i_eng] += (NCG*Int);
  }
  
  ++i_eng;
  
  // adding rho_e_plus_p computed at the cg to each node //

  Int = 0.0;
  for (j=0; j<4; ++j){
    Int += V[nodeNum(j)][0]*((1.0/gam1)*(V[nodeNum(j)][4]/V[nodeNum(j)][0])+0.5*squ[j])+V[nodeNum(j)][4];
  }
  for (j=0; j<4; ++j){
    Eng_Test[nodeNum(j)][i_eng] += (NCG*Int);
  }


  // adding rho_s_dtdxj computed at the cg to each node //

  // step-1: computing temperature at every node
 
  double t[4]; 
  for (j=0; j<4; ++j)
    t[j] = V[nodeNum(j)][4]/(R*V[nodeNum(j)][0]);

  // step-2: computing derivative of temp at the cg 

  double dtdxj[3];  
  dtdxj[0] = (dp1dxj[0][0]*t[0] + dp1dxj[1][0]*t[1] + dp1dxj[2][0]*t[2] + dp1dxj[3][0]*t[3]);
  dtdxj[1] = (dp1dxj[0][1]*t[0] + dp1dxj[1][1]*t[1] + dp1dxj[2][1]*t[2] + dp1dxj[3][1]*t[3]);
  dtdxj[2] = (dp1dxj[0][2]*t[0] + dp1dxj[1][2]*t[1] + dp1dxj[2][2]*t[2] + dp1dxj[3][2]*t[3]);

  // step-3: assembling step
  
  for (i=0; i<3; ++i){
    Int = 0.0;
    for (j=0; j<4; ++j){
      Int += V[nodeNum(j)][0]*sqrt2S2*dtdxj[i];
    }
    ++i_eng;
    for (j=0; j<4; ++j){
      Eng_Test[nodeNum(j)][i_eng] += (NCG*Int);
    }
  }

  // adding dtdxj computed at the cg to each node //
  
  for (i=0; i<3; ++i){
    Int = 0.0;
    for (j=0; j<4; ++j){
      Int += dtdxj[i];
    }
    ++i_eng;  
    for (j=0; j<4; ++j){
      Eng_Test[nodeNum(j)][i_eng] += vol*Int;
    }
  }       

}


//------------------------------------------------------------------------------

template<int dim>
void ElemTet::computeMBarAndM(DynamicVMSTerm *dvmst,
			      SVec<double,dim> **VBar,
			      SVec<double,1> **volRatio,
			      SVec<double,3> &X,
			      SVec<double,dim> &V,
			      SVec<double,dim> &MBar,
			      SVec<double,dim> &M)
{

   int i, j, k;
   const double twothird = 2.0/3.0;
   bool clip = false;

   double dp1dxj[4][3];
   double vol = computeGradientP1Function(X, dp1dxj);

   double *v[4]       = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)], V[nodeNum(3)]};
   double *vbar[4]    = {(*VBar[0])[nodeNum(0)], (*VBar[0])[nodeNum(1)],
                         (*VBar[0])[nodeNum(2)], (*VBar[0])[nodeNum(3)]};
   double *vbarbar[4] = {(*VBar[1])[nodeNum(0)], (*VBar[1])[nodeNum(1)],
                         (*VBar[1])[nodeNum(2)], (*VBar[1])[nodeNum(3)]};

   double vr[4] = {(*volRatio[0])[nodeNum(0)][0], (*volRatio[0])[nodeNum(1)][0],
                   (*volRatio[0])[nodeNum(2)][0], (*volRatio[0])[nodeNum(3)][0]};

   double invVolR[4] = {1.0/vr[0], 1.0/vr[1], 1.0/vr[2], 1.0/vr[3]};

   double r1[3][dim], r2[3][dim];

   double Cs[4] = {1.0,1.0,1.0,1.0};
   double Pt[4] = {1.0,1.0,1.0,1.0};

   dvmst->compute(Cs, Pt, vol, dp1dxj, vbar, v, reinterpret_cast<double *>(r1), X, nodeNum(), clip);

   for (i = 0; i < 4; ++i)
     Cs[i] = pow(invVolR[i], twothird);

   dvmst->compute(Cs, Pt, vol, dp1dxj, vbarbar, vbar, reinterpret_cast<double *>(r2), X, nodeNum(), clip);

   for (int j=0; j<4; ++j) {
     int idx = nodeNum(j);
     for (int k=0; k<dim; ++k) {
       MBar[idx][k] += vol * (r2[0][k] * dp1dxj[j][0]
                            + r2[1][k] * dp1dxj[j][1]
                            + r2[2][k] * dp1dxj[j][2]);
       M[idx][k] += vol * (r1[0][k] * dp1dxj[j][0]
                         + r1[1][k] * dp1dxj[j][1]
                         + r1[2][k] * dp1dxj[j][2]);
     }
   }

}

//-----------------------------------------------------------------------

template<int dim>
void ElemTet::computeDynamicVMSTerm(DynamicVMSTerm *dvmst,
                                SVec<double,dim> **VBar,
                                SVec<double,3> &X,
                                SVec<double,dim> &V,
                                SVec<double,dim> &S,
                                Vec<double> &CsDelSq,
                                Vec<double> &PrT,
                                Vec<double> *Cs,
                                Vec<double> &Delta)

{

  double dp1dxj[4][3];
  double vol = computeGradientP1Function(X, dp1dxj);
  bool   clip = true;     // flag that tiggers clipping of cs and pt values

  double cs[4] = {CsDelSq[nodeNum(0)], CsDelSq[nodeNum(1)], CsDelSq[nodeNum(2)], CsDelSq[nodeNum(3)]};
  double pt[4] = {PrT[nodeNum(0)], PrT[nodeNum(1)], PrT[nodeNum(2)], PrT[nodeNum(3)]};

  double *v[4]    = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)], V[nodeNum(3)]};
  double *vbar[4] = {(*VBar[0])[nodeNum(0)], (*VBar[0])[nodeNum(1)],
                     (*VBar[0])[nodeNum(2)], (*VBar[0])[nodeNum(3)]};

  double r[3][dim];

  dvmst->compute(cs, pt, vol, dp1dxj, vbar, v, reinterpret_cast<double *>(r), X, nodeNum(), clip);

  // reynolds stress flux //

  for (int j=0; j<4; ++j) {
    int idx = nodeNum(j);
    for (int k=0; k<dim; ++k) {
      S[idx][k] += vol *( r[0][k] * dp1dxj[j][0]
                        + r[1][k] * dp1dxj[j][1]
                        + r[2][k] * dp1dxj[j][2]);
    }
  }

  // saving nodal Cs values for post processing //

  if (Cs){
    double Dt[4] =  {Delta[nodeNum(0)], Delta[nodeNum(1)],
                     Delta[nodeNum(2)], Delta[nodeNum(3)]};
    for(int i=0; i<4; ++i) {
       if (cs[i] != 0.0) {
          if(cs[i] < 0.0) (*Cs)[nodeNum(i)] = -sqrt(fabs(cs[i]))/Dt[i];
          else  (*Cs)[nodeNum(i)] = sqrt(cs[i])/Dt[i];
       }
       else  (*Cs)[nodeNum(i)] = 0.0;
    }
  }

}

//-----------------------------------------------------------------------

template<int dim>
void ElemTet::computeVMSLESTerm(VMSLESTerm *vmst,
				SVec<double,dim> &VBar,
				SVec<double,3> &X,
				SVec<double,dim> &V,
				SVec<double,dim> &Sigma)

{

  double dp1dxj[4][3];

  double vol = computeGradientP1Function(X, dp1dxj);

  double *v[4] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)], V[nodeNum(3)]};

  double *vbar[4] = {VBar[nodeNum(0)], VBar[nodeNum(1)], VBar[nodeNum(2)], VBar[nodeNum(3)]};

  double r[3][dim];

  vmst->compute(vol, dp1dxj, vbar, v, reinterpret_cast<double *>(r), X, nodeNum());

  for (int j=0; j<4; ++j) {
    int idx = nodeNum(j);
    for (int k=0; k<dim; ++k) {
      Sigma[idx][k] += vol * ( r[0][k] * dp1dxj[j][0] + r[1][k] * dp1dxj[j][1] +
                               r[2][k] * dp1dxj[j][2] );
    }
  }

}

//-----------------------------------------------------------------------

template<int dim>
void ElemTet::computeSmagorinskyLESTerm(SmagorinskyLESTerm *smag, SVec<double,3> &X,
					SVec<double,dim> &V, SVec<double,dim> &R)
{
  
  double dp1dxj[4][3];
  double vol = computeGradientP1Function(X, dp1dxj);
  double *v[4] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)], V[nodeNum(3)]};
  
  double r[3][dim];
  
  smag->compute(vol, dp1dxj, v, reinterpret_cast<double *>(r), X, nodeNum());
  
  for (int j=0; j<4; ++j) {
    int idx = nodeNum(j);
    for (int k=0; k<dim; ++k) {
      R[idx][k] += vol * ( r[0][k] * dp1dxj[j][0] + r[1][k] * dp1dxj[j][1] +
                           r[2][k] * dp1dxj[j][2] );
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void ElemTet::computeWaleLESTerm(WaleLESTerm *wale, SVec<double,3> &X,
			     SVec<double,dim> &V, SVec<double,dim> &R)
{
  
  double dp1dxj[4][3];
  double vol = computeGradientP1Function(X, dp1dxj);
  double *v[4] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)], V[nodeNum(3)]};
  
  double r[3][dim];
  
  wale->compute(vol, dp1dxj, v, reinterpret_cast<double *>(r), X, nodeNum());
  
  for (int j=0; j<4; ++j) {
    int idx = nodeNum(j);
    for (int k=0; k<dim; ++k) {
      R[idx][k] += vol * ( r[0][k] * dp1dxj[j][0] + r[1][k] * dp1dxj[j][1] +
                           r[2][k] * dp1dxj[j][2] );
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void ElemTet::computeDynamicLESTerm(DynamicLESTerm *dles, SVec<double,2> &Cs,
                                    SVec<double,3> &X, SVec<double,dim> &V, SVec<double,dim> &R)
{

  double dp1dxj[4][3];
  double vol = computeGradientP1Function(X, dp1dxj);
  double *v[4] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)], V[nodeNum(3)]};
  double cs[4] = {Cs[nodeNum(0)][0], Cs[nodeNum(1)][0],
                  Cs[nodeNum(2)][0], Cs[nodeNum(3)][0]};
  double pt[4] = {Cs[nodeNum(0)][1], Cs[nodeNum(1)][1],
                  Cs[nodeNum(2)][1], Cs[nodeNum(3)][1]};

  double r[3][dim];

  dles->compute(vol, dp1dxj, v, cs, pt, reinterpret_cast<double *>(r), X, nodeNum());

  for (int j=0; j<4; ++j) {
    int idx = nodeNum(j);
    for (int k=0; k<dim; ++k) {
      R[idx][k] += vol * ( r[0][k] * dp1dxj[j][0] + r[1][k] * dp1dxj[j][1] +
                           r[2][k] * dp1dxj[j][2] );
    }
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void ElemTet::computeJacobianGalerkinTerm(FemEquationTerm *fet, SVec<double,3> &X, 
					  Vec<double> &ctrlVol, Vec<double> &d2wall, 
					  SVec<double,dim> &V, GenMat<Scalar,neq> &A,
                                          Vec<GhostPoint<dim> *> *ghostPoints, LevelSetStructure *LSS)
{
  // In the case of an embedded simulation, check if the tetrahedra is actually active
  bool isTetInactive=true,isAtTheInterface=true;
  if(ghostPoints) // Then LSS is also a non null pointer. It has already been checked in Domain.C
    {
      for(int i=0;i<4;++i) 
	{
	  isTetInactive    = (isTetInactive && !(LSS->isActive(0,nodeNum(i))));
	  isAtTheInterface = (isAtTheInterface && LSS->isActive(0,nodeNum(i)));
	}
      if(isTetInactive) return;
    }
  isAtTheInterface = !isAtTheInterface;

  // Common Part
  double dp1dxj[4][3];
  double vol = computeGradientP1Function(X, dp1dxj);

  double d2w[4] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)],
                   d2wall[nodeNum(2)], d2wall[nodeNum(3)]};
  double *v[4] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)], V[nodeNum(3)]};

  // Embedded Case. Replace states in the structure by Ghost States
  if(ghostPoints)
  {
   GhostPoint<dim> *gp;
   for(int i=0;i<4;++i)
   {
    gp = ghostPoints->operator[](nodeNum(i)); 
    if(gp) v[i] = gp->getPrimitiveState();
   }
  }

  double dRdU[4][3][neq*neq], dSdU[4][neq*neq], dPdU[4][4][neq*neq];
  bool porousTermExists = fet->computeJacobianVolumeTerm(dp1dxj, d2w, v, reinterpret_cast<double *>(dRdU),
                                                         reinterpret_cast<double *>(dSdU), reinterpret_cast<double *>(dPdU),
                                                         vol, X, nodeNum(), volume_id);

  bool sourceTermExists = fet->doesSourceTermExist();

  dp1dxj[0][0] *= vol;
  dp1dxj[0][1] *= vol;
  dp1dxj[0][2] *= vol;

  dp1dxj[1][0] *= vol;
  dp1dxj[1][1] *= vol;
  dp1dxj[1][2] *= vol;

  dp1dxj[2][0] *= vol;
  dp1dxj[2][1] *= vol;
  dp1dxj[2][2] *= vol;

  dp1dxj[3][0] *= vol;
  dp1dxj[3][1] *= vol;
  dp1dxj[3][2] *= vol;

  double vol4 = vol * fourth;

  // diagonal matrices

  for (int k=0; k<4; ++k) {
    Scalar *Aii = A.getElem_ii(nodeNum(k));
    int m;

    for (m=0; m<neq*neq; ++m)
      Aii[m] += (dRdU[k][0][m] * dp1dxj[k][0] + dRdU[k][1][m] * dp1dxj[k][1] +
                 dRdU[k][2][m] * dp1dxj[k][2]);

    if (sourceTermExists)
      for (m=0; m<neq*neq; ++m)
        Aii[m] -= vol4 * dSdU[k][m];

    if (porousTermExists)
       for (m=0;m<neq*neq;++m)
         Aii[m] += dPdU[k][k][m];

  }

  // off-diagonal matrices

  for (int l=0; l<6; ++l) {
    int i, j;
    if (nodeNum( edgeEnd(l,0) ) < nodeNum( edgeEnd(l,1) )) {
      i = edgeEnd(l,0);
      j = edgeEnd(l,1);
    }
    else {
      i = edgeEnd(l,1);
      j = edgeEnd(l,0);
    }

    if (LSS->isActive(0,nodeNum(i)) && LSS->isActive(0,nodeNum(j))) {
      Scalar *Aij = A.getElem_ij(edgeNum(l));
      Scalar *Aji = A.getElem_ji(edgeNum(l));

      if (Aij && Aji) {
	double cij = 1.0 / ctrlVol[ nodeNum(i) ];
	double cji = 1.0 / ctrlVol[ nodeNum(j) ];
	int m;
	
	for (m=0; m<neq*neq; ++m) {
	  Aij[m] += cij * (dRdU[j][0][m] * dp1dxj[i][0] + dRdU[j][1][m] * dp1dxj[i][1] +
			   dRdU[j][2][m] * dp1dxj[i][2]);
	  Aji[m] += cji * (dRdU[i][0][m] * dp1dxj[j][0] + dRdU[i][1][m] * dp1dxj[j][1] +
			   dRdU[i][2][m] * dp1dxj[j][2]);
	}
	
	if (sourceTermExists) {
	  double cij4 = cij * vol4;
	  double cji4 = cji * vol4;
	  for (m=0; m<neq*neq; ++m) {
	    Aij[m] -= cij4 * dSdU[j][m];
	    Aji[m] -= cji4 * dSdU[i][m];
	  }
	}
	
	if (porousTermExists) {
	  for (m=1;m<neq*neq;++m) {
	    Aij[m] += cij * dPdU[i][j][m];
	    Aji[m] += cji * dPdU[j][i][m];
	  }
	}
      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void ElemTet::computeFaceGalerkinTerm(FemEquationTerm *fet, int face[3], int code, Vec3D &n, 
				      SVec<double,3> &X, Vec<double> &d2wall, double *Vwall, 
				      SVec<double,dim> &V, SVec<double,dim> &R)
{

  double dp1dxj[4][3];
  computeGradientP1Function(X, dp1dxj);

  double d2w[4] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)], 
		   d2wall[nodeNum(2)], d2wall[nodeNum(3)]};
  double *v[4] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)], V[nodeNum(3)]};

  double r[dim];
  fet->computeSurfaceTerm(dp1dxj, code, n, d2w, Vwall, v, r);

  for (int l=0; l<3; ++l)
    for (int k=0; k<dim; ++k)
      R[ face[l] ][k] -= third * r[k];

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void ElemTet::computeDerivativeOfFaceGalerkinTerm(FemEquationTerm *fet, int face[3], int code, Vec3D &n, Vec3D &dn,
				  SVec<double,3> &X, SVec<double,3> &dX, Vec<double> &d2wall, double *Vwall, double *dVwall,
				  SVec<double,dim> &V, SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR)
{

  double dp1dxj[4][3];
  computeGradientP1Function(X, dp1dxj);

  double ddp1dxj[4][3];
  computeDerivativeOfGradientP1Function(X, dX, ddp1dxj);

  double d2w[4] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)],
		   d2wall[nodeNum(2)], d2wall[nodeNum(3)]};
  double *v[4] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)], V[nodeNum(3)]};
  double *dv[4] = {dV[nodeNum(0)], dV[nodeNum(1)], dV[nodeNum(2)], dV[nodeNum(3)]};

  double dr[dim];
  fet->computeDerivativeOfSurfaceTerm(dp1dxj, ddp1dxj, code, n, dn, d2w, Vwall, dVwall, v, dv, dMach, dr);

  for (int l=0; l<3; ++l)
    for (int k=0; k<dim; ++k)
      dR[ face[l] ][k] -= third * dr[k];

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void ElemTet::computeFaceJacobianGalerkinTerm(FemEquationTerm *fet, int face[3], int code, 
					      Vec3D &n, SVec<double,3> &X, Vec<double> &ctrlVol,
					      Vec<double> &d2wall, double *Vwall, 
					      SVec<double,dim> &V, GenMat<Scalar,neq> &A)
{

  double dp1dxj[4][3];
  computeGradientP1Function(X, dp1dxj);

  double d2w[4] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)], 
		   d2wall[nodeNum(2)], d2wall[nodeNum(3)]};
  double *v[4] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)], V[nodeNum(3)]};

  double dRdU[4][neq*neq];
  fet->computeJacobianSurfaceTerm(dp1dxj, code, n, d2w, Vwall, v, 
				  reinterpret_cast<double *>(dRdU));

  for (int k=0; k<4; ++k) {
    if (nodeNum(k) == face[0] || nodeNum(k) == face[1] || nodeNum(k) == face[2]) {
      Scalar *Aii = A.getElem_ii(nodeNum(k));
      for (int m=0; m<neq*neq; ++m)
	Aii[m] -= third * dRdU[k][m];
    }
  }

  for (int l=0; l<6; ++l) {

    int i, j;
    if (nodeNum( edgeEnd(l,0) ) < nodeNum( edgeEnd(l,1) )) {
      i = edgeEnd(l,0);
      j = edgeEnd(l,1);
    } 
    else {
      i = edgeEnd(l,1);
      j = edgeEnd(l,0);
    }

    Scalar *Aij = A.getElem_ij(edgeNum(l));
    Scalar *Aji = A.getElem_ji(edgeNum(l));

    if ( Aij && ( nodeNum(i) == face[0] || 
		  nodeNum(i) == face[1] ||
		  nodeNum(i) == face[2] ) ) {

      double cij = third / ctrlVol[ nodeNum(i) ];
      for (int m=0; m<neq*neq; ++m)
	Aij[m] -= cij * dRdU[j][m];
    }

    if ( Aji && ( nodeNum(j) == face[0] || 
		  nodeNum(j) == face[1] ||
		  nodeNum(j) == face[2] ) ) {

      double cji = third / ctrlVol[ nodeNum(j) ];
      for (int m=0; m<neq*neq; ++m)
	Aji[m] -= cji * dRdU[i][m];
    }

  }

}



//------------------------------------------------------------------------------
//--------------- REINITIALIZATION LEVEL SET    --------------------------------
//------------------------------------------------------------------------------
template<int dim>
int ElemTet::findLSIntersectionPoint(int lsdim, SVec<double,dim> &Phi, SVec<double,dim> &ddx,
                                 SVec<double,dim> &ddy, SVec<double,dim> &ddz,
                                 SVec<double,3> &X,
                                 int reorder[4], Vec3D P[4],
                                 int typeTracking)
{

  // 1 - find which case we are dealing with, ie how many nodes have
  //     positive phis and how many have negative phis
  int positive = 0;
  int negative = 0;
  int zero     = 0;
  for (int i=0; i<4; i++)
    if (Phi[nodeNum(i)][lsdim]<0.0)
      negative++;
    else if(Phi[nodeNum(i)][lsdim]>0.0)
      positive++;
    else
      zero++;
  
  if(negative<1 && positive<1){
    fprintf(stdout, "Error: tetrahedron has only zero phi point values\n");
    exit(1);
  }

  // 2 - orient the tet if necessary (node renumbering from 0 to 3,
  //         which is different again from the local node numbering!)
  //     if all nodes have same phi sign --> nothing to do
  //     if one node is different from the others --> make it be the node 0
  //     if two nodes are different --> first node is unchanged ie reorder[0]=0
  //                                    make sure that reoder[1] has same sign of phi as reorder[0]
	//     if there are nodes with value 0 --> particular cases....

  int scenario;                         // which configuration to run
  int tempi = 0;                        // temporary variables to determine reordering
  int tempi2 = 0;                        // temporary variables to determine reordering

  if((negative==0 || positive==0) && zero==0)
    scenario = 0;
  else if(positive==2 && negative==2){
    scenario = 2;
    if(Phi[nodeNum(0)][lsdim]*Phi[nodeNum(1)][lsdim]<0.0){
    //swap if need be so that reorder[0] and reorder[1] have same sign
      if((Phi[nodeNum(0)][lsdim])*(Phi[nodeNum(2)][lsdim])>0.0){
        reorder[1] = 2;
        reorder[2] = 3;
        reorder[3] = 1;
      }else{
        reorder[1] = 3;
        reorder[2] = 1;
        reorder[3] = 2;
      }
    }
  }
  else if((positive==1 && negative==3) ||
          (positive==3 && negative==1)){//1-vs-3 case not including zero cases
    scenario = 1;
    if(positive==1){ // we want to find i such that Phi[nodeNum(i)][lsdim]>0.0
      while(!(Phi[nodeNum(tempi)][lsdim]>0.0))
        tempi++;
    }
    else if(negative==1){
      while(!(Phi[nodeNum(tempi)][lsdim]<0.0))
        tempi++;
    }

  }
  // cases with zero-valued phis
  else if(zero>0){
    if(zero==1 && (positive==3 || negative==3)){
      scenario = 3;
      while(Phi[nodeNum(tempi)][lsdim]!=0.0)
        tempi++;
    }
    else if(zero==1 && (positive==2 || negative==2)){
      scenario = 4;
      if(positive==1){
        while(!(Phi[nodeNum(tempi)][lsdim]>0.0))
          tempi++;
      }
      if(negative==1){
        while(!(Phi[nodeNum(tempi)][lsdim]<0.0))
          tempi++;
      }
      while(Phi[nodeNum(tempi2)][lsdim]!=0.0)
        tempi2++;
    }
    else if(zero==2 && (positive==2 || negative==2)){
      scenario = 5;
      while(Phi[nodeNum(tempi)][lsdim]!=0.0)
        tempi++;
      tempi2 = tempi+1;
      while(Phi[nodeNum(tempi2)][lsdim]!=0.0)
        tempi2++;
			
    }
    else if(zero==2 && (positive==1 && negative==1)){
      scenario = 6;
      while(!(Phi[nodeNum(tempi)][lsdim]>0.0))
        tempi++;
    }
    else if(zero==3 && (positive==1 || negative==1)){
      scenario = 7;
      while(Phi[nodeNum(tempi)][lsdim]==0.0)
        tempi++;
    }
    else{
      fprintf(stdout, "*** Error: in Tet, # of negative, positive and zero-valued phis are %d %d %d\n", negative, positive, zero);
      exit(1);
    }

  }


	// 3 - reordering of nodes for cases where one node (tempi) was different from other ones
  if(scenario == 1 || scenario == 3 ||
     scenario == 6 || scenario ==7){
    if(tempi==1){
      reorder[0] = 1;
      reorder[1] = 2;
      reorder[2] = 0;
      reorder[3] = 3;
    }else if(tempi==2){
      reorder[0] = 2;
      reorder[1] = 0;
      reorder[2] = 1;
      reorder[3] = 3;
    }else if(tempi==3){
      reorder[0] = 3;
      reorder[1] = 0;
      reorder[2] = 2;
      reorder[3] = 1;
    }
  }
  else if(scenario == 4 ){
    if(tempi==0 && tempi2==1){} //nothing to do
    else if(tempi==0 && tempi2==2){
      reorder[0] = 0;
      reorder[1] = 2;
      reorder[2] = 3;
      reorder[3] = 1;
    }
    else if(tempi==0 && tempi2==3){
      reorder[0] = 0;
      reorder[1] = 3;
      reorder[2] = 1;
      reorder[3] = 2;
    }
    else if(tempi==1 && tempi2==2){
      reorder[0] = 1;
      reorder[1] = 2;
      reorder[2] = 0;
      reorder[3] = 3;
    }
    else if(tempi==1 && tempi2==3){
      reorder[0] = 1;
      reorder[1] = 3;
      reorder[2] = 0;
      reorder[3] = 2;
    }
    else if(tempi==1 && tempi2==0){
      reorder[0] = 1;
      reorder[1] = 0;
      reorder[2] = 3;
      reorder[3] = 2;
    }
    else if(tempi==2 && tempi2==3){
      reorder[0] = 2;
      reorder[1] = 3;
      reorder[2] = 0;
      reorder[3] = 1;
    }
    else if(tempi==2 && tempi2==0){
      reorder[0] = 2;
      reorder[1] = 0;
      reorder[2] = 1;
      reorder[3] = 3;
    }
    else if(tempi==2 && tempi2==1){
      reorder[0] = 2;
      reorder[1] = 1;
      reorder[2] = 3;
      reorder[3] = 0;
    }
    else if(tempi==3 && tempi2==0){
      reorder[0] = 3;
      reorder[1] = 0;
      reorder[2] = 2;
      reorder[3] = 1;
    }
    else if(tempi==3 && tempi2==1){
      reorder[0] = 3;
      reorder[1] = 1;
      reorder[2] = 0;
      reorder[3] = 2;
    }
    else if(tempi==3 && tempi2==2){
      reorder[0] = 3;
      reorder[1] = 2;
      reorder[2] = 1;
      reorder[3] = 0;
    }
  }
  else if(scenario == 5 ){
    if(tempi==0 && tempi2==1){} //nothing to do
    else if(tempi==0 && tempi2==2){
      reorder[0] = 0;
      reorder[1] = 2;
      reorder[2] = 3;
      reorder[3] = 1;
    }
    else if(tempi==0 && tempi2==3){
      reorder[0] = 0;
      reorder[1] = 3;
      reorder[2] = 1;
      reorder[3] = 2;
    }
    else if(tempi==1 && tempi2==2){
      reorder[0] = 1;
      reorder[1] = 2;
      reorder[2] = 0;
      reorder[3] = 3;
    }
    else if(tempi==1 && tempi2==3){
      reorder[0] = 1;
      reorder[1] = 3;
      reorder[2] = 0;
      reorder[3] = 2;
    }
    else if(tempi==2 && tempi2==3){
      reorder[0] = 2;
      reorder[1] = 3;
      reorder[2] = 0;
      reorder[3] = 1;
    }
  }




  // 4 - determination of the coordinates of the intersection points P with 
  //     the material interface
  if(typeTracking == MultiFluidData::LINEAR){
    findLSIntersectionPointLinear(lsdim, Phi,ddx,ddy,ddz,X,reorder,P,scenario);
    return scenario;
  }
  else if(typeTracking == MultiFluidData::GRADIENT){
    findLSIntersectionPointGradient(lsdim,Phi,ddx,ddy,ddz,X,reorder,P,scenario);
    return scenario;
  }
  else if(typeTracking == MultiFluidData::HERMITE){
    int err = findLSIntersectionPointHermite(lsdim,Phi,ddx,ddy,ddz,X,reorder,P,scenario);
    if(err>0) findLSIntersectionPointLinear(lsdim,Phi,ddx,ddy,ddz,X,reorder,P,scenario);
    return scenario;
  }else{
    fprintf(stdout, "Problem in Tet\n");
    exit(1);
  }

}

//------------------------------------------------------------------------------
template<int dim>
void ElemTet::findLSIntersectionPointLinear(int lsdim, SVec<double,dim> &Phi,  SVec<double,dim> &ddx,
                                 SVec<double,dim> &ddy, SVec<double,dim> &ddz,
                                 SVec<double,3> &X,
                                 int reorder[4], Vec3D P[4], int scenario)
{

  Vec3D C0(X[nodeNum(reorder[0])][0],X[nodeNum(reorder[0])][1],X[nodeNum(reorder[0])][2]);
  Vec3D C1(X[nodeNum(reorder[1])][0],X[nodeNum(reorder[1])][1],X[nodeNum(reorder[1])][2]);
  Vec3D C2(X[nodeNum(reorder[2])][0],X[nodeNum(reorder[2])][1],X[nodeNum(reorder[2])][2]);
  Vec3D C3(X[nodeNum(reorder[3])][0],X[nodeNum(reorder[3])][1],X[nodeNum(reorder[3])][2]);

  double ksi[4] = {-1.0, -1.0, -1.0, -1.0};

  // 3 - find the intersection point when they exist
  if (scenario==0){ //sign(phi) is constant in tet
  }
  else if (scenario==2){
  // nodes reorder[0] and reorder[1] have same sign1
  // nodes reorder[2] and reorder[3] have same sign2
  // the plane phi=0 will cross edge reorder[0]-reorder[3] in P3
  //                                 reorder[0]-reorder[2] in P2
  // the plane phi=0 will cross edge reorder[1]-reorder[3] in P1
  //                                 reorder[1]-reorder[2] in P0

    //parametric coordinates of P0, P1, P2, P3
    ksi[0] = Phi[nodeNum(reorder[1])][lsdim]/(Phi[nodeNum(reorder[1])][lsdim]-Phi[nodeNum(reorder[2])][lsdim]);
    ksi[1] = Phi[nodeNum(reorder[1])][lsdim]/(Phi[nodeNum(reorder[1])][lsdim]-Phi[nodeNum(reorder[3])][lsdim]);
    ksi[2] = Phi[nodeNum(reorder[0])][lsdim]/(Phi[nodeNum(reorder[0])][lsdim]-Phi[nodeNum(reorder[2])][lsdim]);
    ksi[3] = Phi[nodeNum(reorder[0])][lsdim]/(Phi[nodeNum(reorder[0])][lsdim]-Phi[nodeNum(reorder[3])][lsdim]);
    P[0] = (1.0-ksi[0])* C1 + ksi[0] * C2;
    P[1] = (1.0-ksi[1])* C1 + ksi[1] * C3;
    P[2] = (1.0-ksi[2])* C0 + ksi[2] * C2;
    P[3] = (1.0-ksi[3])* C0 + ksi[3] * C3;

  }
  else if (scenario==1 || scenario==4 ||
           scenario==6 || scenario==7){
  // node reorder[0] is the only node with sign(phi[reorder[0]]) strictly
  // the plane phi=0 will cross edge reorder[0]-reorder[1] in P0  (maybe in reorder[1])
  // the plane phi=0 will cross edge reorder[0]-reorder[2] in P1  (maybe in reorder[2])
  // the plane phi=0 will cross edge reorder[0]-reorder[3] in P2  (maybe in reorder[3])
  // Note that Pk can be reorder[k] itself (k=0,1,2)

    //parametric coordinates of P0, P1, P2 on their edge
    ksi[0] = Phi[nodeNum(reorder[0])][lsdim]/(Phi[nodeNum(reorder[0])][lsdim]-Phi[nodeNum(reorder[1])][lsdim]);
    ksi[1] = Phi[nodeNum(reorder[0])][lsdim]/(Phi[nodeNum(reorder[0])][lsdim]-Phi[nodeNum(reorder[2])][lsdim]);
    ksi[2] = Phi[nodeNum(reorder[0])][lsdim]/(Phi[nodeNum(reorder[0])][lsdim]-Phi[nodeNum(reorder[3])][lsdim]);

    //physical coordinates of P0, P1, P2
    P[0] = (1.0-ksi[0]) * C0 + ksi[0] * C1;
    P[1] = (1.0-ksi[1]) * C0 + ksi[1] * C2;
    P[2] = (1.0-ksi[2]) * C0 + ksi[2] * C3;
		P[3] = C3;
    //P[3] = C3 is not modified and should not be used later.

  }
  else if (scenario==3){
  // node reorder[0] has phi=0
  // all other nodes have same sign
    P[0] = C0; P[1] = C1; P[2] = C2; P[3] = C3;
  }else if (scenario==5){
  // nodes reorder[0] and reorder[1] have phi=0
  // other nodes have same sign
  //nothing to do for now
    P[0] = C0; P[1] = C1; P[2] = C2; P[3] = C3;
  }

}

//------------------------------------------------------------------------------
template<int dim>
void ElemTet::findLSIntersectionPointGradient(int lsdim, SVec<double,dim> &Phi,  SVec<double,dim> &ddx,
                                 SVec<double,dim> &ddy, SVec<double,dim> &ddz,
                                 SVec<double,3> &X,
                                 int reorder[4], Vec3D P[4], int scenario)
{
  fprintf(stdout, "findLSIntersectionPointGradient\n");
// the variation of phi is not assumed to be linear in the tet.
// we approximate the variations of phi around point i as a linear function,
// for which the zero is found. Same is done for point j. Then the mean of those
// two zeros is considered as the intersection of the phi=0 plane and the edges
// considered, ie i-j

  Vec3D C0(X[nodeNum(reorder[0])][0],X[nodeNum(reorder[0])][1],X[nodeNum(reorder[0])][2]);
  Vec3D C1(X[nodeNum(reorder[1])][0],X[nodeNum(reorder[1])][1],X[nodeNum(reorder[1])][2]);
  Vec3D C2(X[nodeNum(reorder[2])][0],X[nodeNum(reorder[2])][1],X[nodeNum(reorder[2])][2]);
  Vec3D C3(X[nodeNum(reorder[3])][0],X[nodeNum(reorder[3])][1],X[nodeNum(reorder[3])][2]);
  Vec3D C[4] = {C0,C1,C2,C3};

  if(scenario==0){
  // nothing to do since sign(phi) = constant in tet

  }else if(scenario==1){
  // node reorder[0] is the only node with sign(phi[reorder[0]]) strictly
  // the plane phi=0 will cross edge reorder[0]-reorder[1] in P0
  // the plane phi=0 will cross edge reorder[0]-reorder[2] in P1
  // the plane phi=0 will cross edge reorder[0]-reorder[3] in P2
  // Note that Pk can be reorder[k] itself (k=0,1,2)
    double phii,phij,gradi,gradj;
    Vec3D nedge,dphij;
    Vec3D dphii(ddx[nodeNum(reorder[0])][lsdim],ddy[nodeNum(reorder[0])][lsdim],ddz[nodeNum(reorder[0])][lsdim]);

    for(int j=0; j<3; j++){
      nedge = C[j+1] - C[0];
      dphij[0] = ddx[nodeNum(reorder[j+1])][lsdim];
      dphij[1] = ddy[nodeNum(reorder[j+1])][lsdim];
      dphij[2] = ddz[nodeNum(reorder[j+1])][lsdim];
      gradi = dphii*nedge;
      gradj = dphij*nedge;
      phii = Phi[nodeNum(reorder[  0])][lsdim]/gradi;
      phij = Phi[nodeNum(reorder[j+1])][lsdim]/gradj;
      P[j] = 0.5*(C[0] + C[j+1] - (phii+phij)*nedge);
    }


  }else if(scenario==2){
  // nodes reorder[0] and reorder[1] have same sign1
  // nodes reorder[2] and reorder[3] have same sign2
  // the plane phi=0 will cross edge reorder[0]-reorder[3] in P3
  //                                 reorder[0]-reorder[2] in P2
  // the plane phi=0 will cross edge reorder[1]-reorder[3] in P1
  //                                 reorder[1]-reorder[2] in P0
    double phii,phij,gradi,gradj;
    Vec3D nedge,dphij,dphii;

    for(int j=2; j<4; j++){
      nedge = C[j] - C[0];
      dphii[0] = ddx[nodeNum(reorder[0])][lsdim];
      dphii[1] = ddy[nodeNum(reorder[0])][lsdim];
      dphii[2] = ddz[nodeNum(reorder[0])][lsdim];
      dphij[0] = ddx[nodeNum(reorder[j])][lsdim];
      dphij[1] = ddy[nodeNum(reorder[j])][lsdim];
      dphij[2] = ddz[nodeNum(reorder[j])][lsdim];
      gradi = dphii*nedge;
      gradj = dphij*nedge;
      phii = Phi[nodeNum(reorder[  0])][lsdim]/gradi;
      phij = Phi[nodeNum(reorder[  j])][lsdim]/gradj;
      P[j] = 0.5*(C[0] + C[j] - (phii+phij)*nedge);
    }
    for(int j=2; j<4; j++){
      nedge = C[j] - C[1];
      dphii[0] = ddx[nodeNum(reorder[1])][lsdim];
      dphii[1] = ddy[nodeNum(reorder[1])][lsdim];
      dphii[2] = ddz[nodeNum(reorder[1])][lsdim];
      dphij[0] = ddx[nodeNum(reorder[j])][lsdim];
      dphij[1] = ddy[nodeNum(reorder[j])][lsdim];
      dphij[2] = ddz[nodeNum(reorder[j])][lsdim];
      gradi = dphii*nedge;
      gradj = dphij*nedge;
      phii = Phi[nodeNum(reorder[  1])][lsdim]/gradi;
      phij = Phi[nodeNum(reorder[  j])][lsdim]/gradj;
      P[j-2] = 0.5*(C[1] + C[j] - (phii+phij)*nedge);
    }

  }

}
//------------------------------------------------------------------------------
template<int dim>
int ElemTet::findLSIntersectionPointHermite(int lsdim, SVec<double,dim> &Phi,  SVec<double,dim> &ddx,
                                 SVec<double,dim> &ddy, SVec<double,dim> &ddz,
                                 SVec<double,3> &X,
                                 int reorder[4], Vec3D P[4], int scenario)
{
/* the variations of phi are not known in the tet, but we know the values
** of phi at each node as well as the derivatives.
** Thus we use Hermite interpolations with Hermite polynomials of degree 3.
** And on each edge, the root of the Hermite polynomial is found.
** This root is the location of the interface.
** 
** Along each edge that the interface crosses, we consider a function f
** that takes values f1 = f(x1) = phi(x1) and f2 = f(x2) = phi(x2)
** and that has derivatives fp1 = f'(x1) = grad(phi(x1)).unitary(x2x1)
**                          fp2 = f'(x2) = grad(phi(x2)).unitary(x2x1)
** This determines a unique polynomial of degree 3 or less for which
** we need to find a real root between x1 and x2.
** IF there are more than one of these roots, an error is returned and
** interface location is switched back to linear interpolation.
*/

  int count = 0;
  Vec3D C0(X[nodeNum(reorder[0])][0],X[nodeNum(reorder[0])][1],X[nodeNum(reorder[0])][2]);
  Vec3D C1(X[nodeNum(reorder[1])][0],X[nodeNum(reorder[1])][1],X[nodeNum(reorder[1])][2]);
  Vec3D C2(X[nodeNum(reorder[2])][0],X[nodeNum(reorder[2])][1],X[nodeNum(reorder[2])][2]);
  Vec3D C3(X[nodeNum(reorder[3])][0],X[nodeNum(reorder[3])][1],X[nodeNum(reorder[3])][2]);
  Vec3D C[4] = {C0,C1,C2,C3};

  double ksi=-1.0;

  // 3 - find the intersection point when they exist
  if (scenario==0){ //sign(phi) is constant in tet
  }else if(scenario==1){
  // node reorder[0] is the only node with sign(phi[reorder[0]]) strictly
  // the plane phi=0 will cross edge reorder[0]-reorder[1] in P0
  // the plane phi=0 will cross edge reorder[0]-reorder[2] in P1
  // the plane phi=0 will cross edge reorder[0]-reorder[3] in P2
  // Note that Pk can be reorder[k] itself (k=0,1,2)
    double f1,f2,fp1,fp2;
    Vec3D nedge,dphij;
    Vec3D dphii(ddx[nodeNum(reorder[0])][lsdim],ddy[nodeNum(reorder[0])][lsdim],ddz[nodeNum(reorder[0])][lsdim]);

    for(int j=0; j<3; j++){
      nedge = C[j+1] - C[0];
      dphij[0] = ddx[nodeNum(reorder[j+1])][lsdim];
      dphij[1] = ddy[nodeNum(reorder[j+1])][lsdim];
      dphij[2] = ddz[nodeNum(reorder[j+1])][lsdim];
      fp1 = dphii*nedge;
      fp2 = dphij*nedge;
      f1 = Phi[nodeNum(reorder[  0])][lsdim];
      f2 = Phi[nodeNum(reorder[j+1])][lsdim];
      //ksi = findRootPolynomialNewtonRaphson(f1,f2,fp1,fp2);
      count = findRootPolynomialLaguerre(f1,f2,fp1,fp2,ksi);
      if(count>1) return 1;
      P[j] = C[0] + ksi*nedge;
    }


  }else if (scenario==2){
  // nodes reorder[0] and reorder[1] have same sign1
  // nodes reorder[2] and reorder[3] have same sign2
  // the plane phi=0 will cross edge reorder[0]-reorder[3] in P3
  //                                 reorder[0]-reorder[2] in P2
  // the plane phi=0 will cross edge reorder[1]-reorder[3] in P1
  //                                 reorder[1]-reorder[2] in P0
    double f1,f2,fp1,fp2;
    Vec3D nedge,dphij,dphii;

    dphii[0] = ddx[nodeNum(reorder[0])][lsdim];
    dphii[1] = ddy[nodeNum(reorder[0])][lsdim];
    dphii[2] = ddz[nodeNum(reorder[0])][lsdim];
    for(int j=2; j<4; j++){
      nedge = C[j] - C[0];
      dphij[0] = ddx[nodeNum(reorder[j])][lsdim];
      dphij[1] = ddy[nodeNum(reorder[j])][lsdim];
      dphij[2] = ddz[nodeNum(reorder[j])][lsdim];
      fp1 = dphii*nedge;
      fp2 = dphij*nedge;
      f1 = Phi[nodeNum(reorder[  0])][lsdim];
      f2 = Phi[nodeNum(reorder[  j])][lsdim];
      //ksi = findRootPolynomialNewtonRaphson(f1,f2,fp1,fp2);
      count = findRootPolynomialLaguerre(f1,f2,fp1,fp2,ksi);
      if(count>1) return 1;
      P[j] = C[0] + ksi*nedge;
    }
    dphii[0] = ddx[nodeNum(reorder[1])][lsdim];
    dphii[1] = ddy[nodeNum(reorder[1])][lsdim];
    dphii[2] = ddz[nodeNum(reorder[1])][lsdim];
    for(int j=2; j<4; j++){
      nedge = C[j] - C[1];
      dphij[0] = ddx[nodeNum(reorder[j])][lsdim];
      dphij[1] = ddy[nodeNum(reorder[j])][lsdim];
      dphij[2] = ddz[nodeNum(reorder[j])][lsdim];
      fp1 = dphii*nedge;
      fp2 = dphij*nedge;
      f1 = Phi[nodeNum(reorder[  1])][lsdim];
      f2 = Phi[nodeNum(reorder[  j])][lsdim];
      //ksi = findRootPolynomialNewtonRaphson(f1,f2,fp1,fp2);
      count = findRootPolynomialLaguerre(f1,f2,fp1,fp2,ksi);
      if(count>1) return 1;
      P[j-2] = C[1] + ksi*nedge;
    }

  }
  return 0;

}
//------------------------------------------------------------------------------
// PDE resolution for reinitialization of Level set

template<int dim>
void ElemTet::computeDistanceToInterface(int type, SVec<double,3> &X, int reorder[4],
                                Vec3D P[4], SVec<double,dim> &Psi, Vec<int> &Tag)
{
  double psi = -1.0;
  int tag = 0;
  double phi[3] = {0.0,0.0,0.0};
  Vec3D Y0,Y1,Y2; //initialized to zero-vector

  //computation of distance to interface phi=0 in the tet for each point
  //each scenario is treated separately and the sign of Tag (which is
  //necessarily 1) is changed if the new psi value (ie Psi is updated
  //by a lower value psi) is attained at a boundary (edge or vertex)
  //so that it can be treated later correctly



  if(type==1){  //3 intersection points are on edges strictly
    Y1 = P[1]-P[0];
    Y2 = P[2]-P[0];
    for(int i=0; i<4; i++){
      Vec3D XX(X[nodeNum(i)]);
      Y0 = XX-P[0];
      tag = computeDistanceToAll(phi,Y0,Y1,Y2,psi);
      if(psi<Psi[nodeNum(i)][0]){
        Psi[nodeNum(i)][0] = psi;
        if(Tag[nodeNum(i)]==-1) Tag[nodeNum(i)]=tag;
      }
    }
  }	
  if(type==2){  //4 intersection points are on edges strictly
    Y1 = P[1]-P[0];
    Y2 = P[2]-P[0];
    for (int i=0; i<4; i++){
      Vec3D XX(X[nodeNum(i)]);
      Y0 = XX-P[0];
      tag = computeDistanceToAll(phi,Y0,Y1,Y2,psi);
      if(psi<Psi[nodeNum(i)][0]){
        Psi[nodeNum(i)][0] = psi;
        if(Tag[nodeNum(i)]==-1) Tag[nodeNum(i)]=tag;
      }
    }
    Y1 = P[1]-P[3];
    Y2 = P[2]-P[3];
    for (int i=0; i<4; i++){
      Vec3D XX(X[nodeNum(i)]);
      Y0 = XX-P[3];
      tag = computeDistanceToAll(phi,Y0,Y1,Y2,psi);
	if(psi<Psi[nodeNum(i)][0]){
        Psi[nodeNum(i)][0] = psi;
        if(Tag[nodeNum(i)]==-1) Tag[nodeNum(i)]=tag;
      }
    }
  }	
  if(type==3){  //1 zero-node and 3 nodes of same sign
    Psi[nodeNum(reorder[0])][0] = 0.0;
    Tag[nodeNum(reorder[0])] = 1;
    for (int i=1; i<4; i++){
      Y0 = P[i] - P[0];
      computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi);
      Psi[nodeNum(reorder[i])][0] = min(psi,Psi[nodeNum(reorder[i])][0]);
    }
		
  }	
  if(type==4){  //1 zero-node and 2 nodes of same sign 
    Y1 = P[1]-P[0];
    Y2 = P[2]-P[0];
    for(int i=0; i<4; i++){
      Vec3D XX(X[nodeNum(i)]);
      Y0 = XX-P[0];
      tag = computeDistanceToAll(phi,Y0,Y1,Y2,psi);
      if(psi<Psi[nodeNum(i)][0]){
        Psi[nodeNum(i)][0] = psi;
        if(Tag[nodeNum(i)]==-1) Tag[nodeNum(i)]=tag;
      }
    }
    Psi[nodeNum(reorder[1])][0] = 0.0;
    Tag[nodeNum(reorder[1])]    = 1;
  }	
  if(type==5){  //2 zero-node and 2 nodes of same sign 
    Psi[nodeNum(reorder[0])][0] = 0.0;
    Tag[nodeNum(reorder[0])]    = 1;
    Psi[nodeNum(reorder[1])][0] = 0.0;
    Tag[nodeNum(reorder[1])]    = 1;

    Y1 = P[1]-P[0];
    for(int i=2; i<4; i++){
      Y0 = P[i]-P[0];
      //psi is overwritten in this function, so that node reorder[3] won t get
      //value from node reorder[2]
      computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi);
      computeDistancePlusPhiToEdge(phi[0],phi[1],Y0,Y1,psi);
      Psi[nodeNum(reorder[i])][0] = min(psi,Psi[nodeNum(reorder[i])][0]);
    }

  }	
  if(type==6){  //2 zero-node and 2 nodes of different signs
    Y1 = P[1]-P[0];
    Y2 = P[2]-P[0];
    for(int i=0; i<4; i++){
      Vec3D XX(X[nodeNum(i)]);
      Y0 = XX-P[0];
      tag = computeDistanceToAll(phi,Y0,Y1,Y2,psi);
      if(psi<Psi[nodeNum(i)][0]){
        Psi[nodeNum(i)][0] = psi;
        if(Tag[nodeNum(i)]==-1) Tag[nodeNum(i)]=tag;
      }
    }
  }	
  if(type==7){  //3 zero-node
    Psi[nodeNum(reorder[1])][0] = 0.0;
    Tag[nodeNum(reorder[1])]    = 1;
    Psi[nodeNum(reorder[2])][0] = 0.0;
    Tag[nodeNum(reorder[2])]    = 1;
    Psi[nodeNum(reorder[3])][0] = 0.0;
    Tag[nodeNum(reorder[3])]    = 1;
    Y1 = P[1]-P[0];
    Y2 = P[2]-P[0];
    Vec3D XX(X[nodeNum(reorder[0])]);
    Y0 = XX-P[0];
    tag = computeDistanceToAll(phi,Y0,Y1,Y2,psi);
    if(psi<Psi[nodeNum(reorder[0])][0]){
      Psi[nodeNum(reorder[0])][0] = psi;
      if(Tag[nodeNum(reorder[0])]==-1) Tag[nodeNum(reorder[0])]=tag;
    }
  }	

}
//------------------------------------------------------------------------------
template<int dim>
void ElemTet::recomputeDistanceToInterface(int type, SVec<double,3> &X, int reorder[4],
                                Vec3D P[4], SVec<double,dim> &Psi, Vec<int> &Tag)
{
  bool found = false;
  double psi = -1.0;
  int tag = 0;
  double phi[3] = {0.0,0.0,0.0};
  Vec3D Y0,Y1,Y2; //initialized to zero-vector
  // we recompute distances differently for nodes that have Tag = -1
  // only when the intersection type is 1, 2, 3, 4 or 5, can that distance
  // be recomputed and only for certain nodes of the tet. 
  // In other cases, there would be no improvement.
  // To recompute that distance, we do not compute the distance from
  // the point to the interface phi=0, but min(phi(x)+|Ax|) where A
  // is the node we consider and x is a point on the chunk of surface
  // opposite to A and with the same sign of phi as A.
  bool show = false;

  if(type==1){  //3 intersection points are on edges strictly
    Vec3D C1(X[nodeNum(reorder[1])]);
    Vec3D C2(X[nodeNum(reorder[2])]);
    Vec3D C3(X[nodeNum(reorder[3])]);
    //node reorder[1]-oppface=P[2]C2C3+P[1]C2P[2]
    if(Tag[nodeNum(reorder[1])]==-1){
    phi[0] = 0.0;
    phi[1] = Psi[nodeNum(reorder[2])][0];
    phi[2] = Psi[nodeNum(reorder[3])][0];
    Y0 = C1-P[2];
    Y1 = C2-P[2];
    Y2 = C3-P[2];
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[1])][0]){
      Psi[nodeNum(reorder[1])][0] = psi;
      if(found) Tag[nodeNum(reorder[1])]    = 1;
    }

    phi[0] = 0.0;
    phi[1] = Psi[nodeNum(reorder[2])][0];
    phi[2] = 0.0;
    Y0 = C1-P[1];
    Y1 = C2-P[1];
    Y2 = P[2]-P[1];
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[1])][0]){
      Psi[nodeNum(reorder[1])][0] = psi;
      if(found) Tag[nodeNum(reorder[1])]    = 1;
    }
    }

    //node reorder[2]-oppface=P[0]C3C1+P[2]C3P[0]
    if(Tag[nodeNum(reorder[2])]==-1){
    phi[0] = 0.0;
    phi[1] = Psi[nodeNum(reorder[3])][0];
    phi[2] = Psi[nodeNum(reorder[1])][0];
    Y0 = C2-P[0];
    Y1 = C3-P[0];
    Y2 = C1-P[0];
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[2])][0]){
      Psi[nodeNum(reorder[2])][0] = psi;
      if(found) Tag[nodeNum(reorder[2])]    = 1;
    }

    phi[0] = 0.0;
    phi[1] = Psi[nodeNum(reorder[3])][0];
    phi[2] = 0.0;
    Y0 = C2-P[2];
    Y1 = C3-P[2];
    Y2 = P[0]-P[2];
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[2])][0]){
      Psi[nodeNum(reorder[2])][0] = psi;
      if(found) Tag[nodeNum(reorder[2])]    = 1;
    }
    }
    //node reorder[3]-oppface=P[1]C1C2+P[0]C1P[1]
    if(Tag[nodeNum(reorder[3])]==-1){
    phi[0] = 0.0;
    phi[1] = Psi[nodeNum(reorder[1])][0];
    phi[2] = Psi[nodeNum(reorder[2])][0];
    Y0 = C3-P[1];
    Y1 = C1-P[1];
    Y2 = C2-P[1];
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[3])][0]){
      Psi[nodeNum(reorder[3])][0] = psi;
      if(found) Tag[nodeNum(reorder[3])]    = 1;
    }

    phi[0] = 0.0;
    phi[1] = Psi[nodeNum(reorder[1])][0];
    phi[2] = 0.0;
    Y0 = C3-P[0];
    Y1 = C1-P[0];
    Y2 = P[1]-P[0];
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[3])][0]){
      Psi[nodeNum(reorder[3])][0] = psi;
      if(found) Tag[nodeNum(reorder[3])]    = 1;
    }
    }
  }
  if(type==2){  //4 intersection points are on edges strictly
    Vec3D C0(X[nodeNum(reorder[0])]);
    Vec3D C1(X[nodeNum(reorder[1])]);
    Vec3D C2(X[nodeNum(reorder[2])]);
    Vec3D C3(X[nodeNum(reorder[3])]);
    //node reorder[0]-oppface=C1P[1]P[0]
    if(Tag[nodeNum(reorder[0])]==-1){
    phi[0] = Psi[nodeNum(reorder[1])][0];
    phi[1] = -phi[0];
    phi[2] = -phi[0];
    Y0 = C0-C1;
    Y1 = P[1]-C1;
    Y2 = P[0]-C1;
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[0])][0]){
      Psi[nodeNum(reorder[0])][0] = psi;
      if(found) Tag[nodeNum(reorder[0])]    = 1;
    }
    }
    //node reorder[1]-oppface=C0P[2]P[3]
    if(Tag[nodeNum(reorder[1])]==-1){
    phi[0] = Psi[nodeNum(reorder[0])][0];
    phi[1] = -phi[0];
    phi[2] = -phi[0];
    Y0 = C1-C0;
    Y1 = P[2]-C0;
    Y2 = P[3]-C0;
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[1])][0]){
      Psi[nodeNum(reorder[1])][0] = psi;
      if(found) Tag[nodeNum(reorder[1])]    = 1;
    }
    }
    //node reorder[2]-oppface=C3P[1]P[3]
    if(Tag[nodeNum(reorder[2])]==-1){
    phi[0] = Psi[nodeNum(reorder[3])][0];
    phi[1] = -phi[0];
    phi[2] = -phi[0];
    Y0 = C2-C3;
    Y1 = P[1]-C3;
    Y2 = P[3]-C3;
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[2])][0]){
      Psi[nodeNum(reorder[2])][0] = psi;
      if(found) Tag[nodeNum(reorder[2])]    = 1;
    }
    }
    //node reorder[3]-oppface=C2P[2]P[0]
    if(Tag[nodeNum(reorder[3])]==-1){
    phi[0] = Psi[nodeNum(reorder[2])][0];
    phi[1] = -phi[0];
    phi[2] = -phi[0];
    Y0 = C3-C2;
    Y1 = P[2]-C2;
    Y2 = P[0]-C2;
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[3])][0]){
      Psi[nodeNum(reorder[3])][0] = psi;
      if(found) Tag[nodeNum(reorder[3])]    = 1;
    }
    }
  }
  if(type==3){  //1 zero-node and 3 nodes of same sign
    //node reorder[1]-oppface=P[0]P[2]P[3]
    if(Tag[nodeNum(reorder[1])]==-1){
    phi[0] = 0.0;
    phi[1] = Psi[nodeNum(reorder[2])][0];
    phi[2] = Psi[nodeNum(reorder[3])][0];
    Y0 = P[1]-P[0];
    Y1 = P[2]-P[0];
    Y2 = P[3]-P[0];
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[1])][0]){
      Psi[nodeNum(reorder[1])][0] = psi;
      if(found) Tag[nodeNum(reorder[1])]    = 1;
    }
    }
    //node reorder[2]-oppface=P[0]P[3]P[1]
    if(Tag[nodeNum(reorder[2])]==-1){
    phi[0] = 0.0;
    phi[1] = Psi[nodeNum(reorder[3])][0];
    phi[2] = Psi[nodeNum(reorder[1])][0];
    Y0 = P[2]-P[0];
    Y1 = P[3]-P[0];
    Y2 = P[1]-P[0];
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[2])][0]){
      Psi[nodeNum(reorder[2])][0] = psi;
      if(found) Tag[nodeNum(reorder[2])]    = 1;
    }
    }
    //node reorder[3]-oppface=P[0]P[1]P[2]
    if(Tag[nodeNum(reorder[3])]==-1){
    phi[0] = 0.0;
    phi[1] = Psi[nodeNum(reorder[1])][0];
    phi[2] = Psi[nodeNum(reorder[2])][0];
    Y0 = P[3]-P[0];
    Y1 = P[1]-P[0];
    Y2 = P[2]-P[0];
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[3])][0]){
      Psi[nodeNum(reorder[3])][0] = psi;
      if(found) Tag[nodeNum(reorder[3])]    = 1;
    }
    }
  }
  if(type==4){  //1 zero-node and 2 nodes of same sign 
    Vec3D C2(X[nodeNum(reorder[2])]);
    Vec3D C3(X[nodeNum(reorder[3])]);
    //node reorder[2]-oppface=P[0]P[2]C3
    if(Tag[nodeNum(reorder[2])]==-1){
    phi[0] = 0.0;
    phi[1] = 0.0;
    phi[2] = Psi[nodeNum(reorder[3])][0];
    Y0 = C2-P[0];
    Y1 = P[2]-P[0];
    Y2 = C3-P[0];
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[2])][0]){
      Psi[nodeNum(reorder[2])][0] = psi;
      if(found) Tag[nodeNum(reorder[2])]    = 1;
    }
    }
    //node reorder[3]-oppface=P[0]C2 P[1]
    if(Tag[nodeNum(reorder[3])]==-1){
    phi[0] = 0.0;
    phi[1] = Psi[nodeNum(reorder[2])][0];
    phi[2] = 0.0;
    Y0 = C3-P[0];
    Y1 = C2-P[0];
    Y2 = P[1]-P[0];
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[3])][0]){
      Psi[nodeNum(reorder[3])][0] = psi;
      if(found) Tag[nodeNum(reorder[3])]    = 1;
    }
    }
  }
  if(type==5){  //2 zero-node and 2 nodes of same sign 
    //node reorder[2]-oppface=P[0]P[3]P[1]
    if(Tag[nodeNum(reorder[2])]==-1){
    phi[0] = Psi[nodeNum(reorder[0])][0];
    phi[1] = Psi[nodeNum(reorder[3])][0]-phi[0];
    phi[2] = Psi[nodeNum(reorder[1])][0]-phi[0];
    Y0 = P[2]-P[0];
    Y1 = P[3]-P[0];
    Y2 = P[1]-P[0];
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[2])][0]){
      Psi[nodeNum(reorder[2])][0] = psi;
      if(found) Tag[nodeNum(reorder[2])]    = 1;
    }
    }

    //node reorder[3]-oppface=P[0]P[1]P[2]
    if(Tag[nodeNum(reorder[3])]==-1){
    phi[0] = Psi[nodeNum(reorder[0])][0];
    phi[1] = Psi[nodeNum(reorder[1])][0]-phi[0];
    phi[2] = Psi[nodeNum(reorder[2])][0]-phi[0];
    Y0 = P[3]-P[0];
    Y1 = P[1]-P[0];
    Y2 = P[2]-P[0];
    computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi,show);
    found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi,show);
    computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi,show);
    if(psi<Psi[nodeNum(reorder[3])][0]){
      Psi[nodeNum(reorder[3])][0] = psi;
      if(found) Tag[nodeNum(reorder[3])]    = 1;
    }
    }
  }
}
//------------------------------------------------------------------------------
template<int dim>
void ElemTet::computeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                    SVec<double,dim> &ddx, SVec<double,dim> &ddy,
                                    SVec<double,dim> &ddz,
                                    SVec<double,dim> &Phi,SVec<double,1> &Psi)
{
  if (!(fabs(Tag[nodeNumTet[0]])==1 && fabs(Tag[nodeNumTet[1]])==1 &&
      fabs(Tag[nodeNumTet[2]])==1 && fabs(Tag[nodeNumTet[3]])==1))
    return;

  // We want to get the values of Psi for the nodes that are closest to 
  // the interface, ie those tagged with value 1

  //find what kind of tetrahedron this is:
  //  0 - no levelset phi=0 in it
  //  1 - levelset phi=0 separates tet in 1 and 3 nodes
  //  2 - levelset phi=0 separates tet in 2 and 2 nodes
  int reorder[4] = {0,1,2,3}; //no change in ordering
  Vec3D P[4] = {X[nodeNum(0)],X[nodeNum(1)],X[nodeNum(2)],X[nodeNum(3)]};

  int type = findLSIntersectionPoint(lsdim, Phi,ddx,ddy,ddz,X,reorder,P,MultiFluidData::LINEAR);
  if(type>0)  computeDistanceToInterface(type,X,reorder,P,Psi,Tag);
}
//------------------------------------------------------------------------------
template<int dim>
void ElemTet::recomputeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                    SVec<double,dim> &ddx, SVec<double,dim> &ddy,
                                    SVec<double,dim> &ddz,
                                    SVec<double,dim> &Phi,SVec<double,1> &Psi)
{
  if (!(Tag[nodeNumTet[0]]==-1 || Tag[nodeNumTet[1]]==-1 ||
      Tag[nodeNumTet[2]]==-1 || Tag[nodeNumTet[3]]==-1))
    return;

  int reorder[4] = {0,1,2,3}; //no change in ordering
  Vec3D P[4] = {X[nodeNum(0)],X[nodeNum(1)],X[nodeNum(2)],X[nodeNum(3)]};

  int type = findLSIntersectionPoint(lsdim,Phi,ddx,ddy,ddz,X,reorder,P,MultiFluidData::LINEAR);
  if(type>0) recomputeDistanceToInterface(type,X,reorder,P,Psi,Tag);
}
//------------------------------------------------------------------------------
template<int dim>
void ElemTet::computeDistanceLevelNodes(int lsdim, Vec<int> &Tag, int level,
                                    SVec<double,3> &X, SVec<double,1> &Psi, SVec<double,dim> &Phi)
{
  if (!(Tag[nodeNumTet[0]]==level || Tag[nodeNumTet[1]]==level ||
       Tag[nodeNumTet[2]]==level || Tag[nodeNumTet[3]]==level   ) )
    return;

  // We want to get the values of Psi for the nodes that are tagged 
  // with values 'level', where level>1 
  // cf Mut, Buscaglia
  for (int i=0; i<4; i++){
    if(Tag[nodeNum(i)]==level||Tag[nodeNum(i)]==level+1){ // compute new value
      double psi = computeDistancePlusPhi(i,X,Psi);
      Psi[nodeNum(i)][0] = min(Psi[nodeNum(i)][0], psi);
    }
  }
}

//------------------------------------------------------------------------------
template<int dim>
double ElemTet::computeDistancePlusPhi(int i, SVec<double,3> &X, SVec<double,dim> &Psi)
{
  // this function computes the following function
  // min(phi(x)+dist(X[nodeNum(i)] - x), x in opposing face to nodeNum(i))
  // to compute this minimum, the location of zero gradient is found inside
  // the tet. If not found, we look at the boundaries, first edges, then 
  // vertices.

  bool show = false;

  double minimum = -1.0; // this function will be a distance eventually (>0.0)
  bool   found = false;

  // list of nodes that define the opposite face.
  int oppi = 3 - i;
  int oppn[3]  = {faceDefTet[oppi][0],faceDefTet[oppi][1],faceDefTet[oppi][2]};
  oppn[0]=nodeNum(oppn[0]); oppn[1]=nodeNum(oppn[1]); oppn[2]=nodeNum(oppn[2]);

  // setup for computations
  double phi[3] = {Psi[oppn[0]][0],
                   Psi[oppn[1]][0]-Psi[oppn[0]][0],
                   Psi[oppn[2]][0]-Psi[oppn[0]][0]};
  // coordinate basis to do our computations (not orthogonal!!)
  Vec3D Y0(X[nodeNum(i)][0]-X[oppn[0]][0],X[nodeNum(i)][1]-X[oppn[0]][1], X[nodeNum(i)][2]-X[oppn[0]][2]);
  Vec3D Y1(X[oppn[1]][0]-X[oppn[0]][0],X[oppn[1]][1]-X[oppn[0]][1], X[oppn[1]][2]-X[oppn[0]][2]);
  Vec3D Y2(X[oppn[2]][0]-X[oppn[0]][0],X[oppn[2]][1]-X[oppn[0]][1], X[oppn[2]][2]-X[oppn[0]][2]);

  computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,minimum,show);

  found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,minimum,show);
	
  computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,minimum,show);

  return minimum;

}

template<int dim, class Obj>
void ElemTet::integrateFunction(Obj* obj,SVec<double,3> &X,SVec<double,dim>& V, void (Obj::*F)(int node, const double* loc,double* f),
				int npt) {

  double vol = computeVolume(X);
  const double locs[5][5] = { {0.0, 0.0,0.0,0.0,0.0},
			      {0.5773502691896257645091488, -0.5773502691896257645091488, 0.0, 0.0, 0.0},
			      {0.7745966692414833770358531, 0.0, -0.7745966692414833770358531, 0.0, 0.0},
			      {0.8611363115940525752239465,0.3399810435848562648026658,-0.3399810435848562648026658,-0.8611363115940525752239465,0.0},
			      {0.9061798459386639927976269,0.5384693101056830910363144,0.0,-0.5384693101056830910363144,-0.9061798459386639927976269} };
  
  const double wgts[5][5] = { { 2.0, 0.0,0.0,0.0,0.0},
			      {1.0,1.0,0.0,0.0,0.0},
			      {0.5555555555555555555555556,0.8888888888888888888888889,0.5555555555555555555555556,0.0,0.0},
			      {0.3478548451374538573730639,0.6521451548625461426269361,0.6521451548625461426269361,0.3478548451374538573730639,0.0},
			      {0.2369268850561890875142640,0.4786286704993664680412915,0.5688888888888888888888889,0.4786286704993664680412915,0.2369268850561890875142640} };

  // First loop through the nodes of the tet
  Vec3D centroid;
  for (int k = 0; k < 3; ++k)
    centroid[k] = 0.25*(X[nodeNum(0)][k]+X[nodeNum(1)][k]+X[nodeNum(2)][k]+X[nodeNum(3)][k]);
  Vec3D faceCnt[4];
  Vec3D edgeCnt[6];
  
  for (int k = 0; k < 3; ++k) {
    edgeCnt[0][k] = 0.5*(X[nodeNum(0)][k]+X[nodeNum(1)][k]);
    edgeCnt[1][k] = 0.5*(X[nodeNum(0)][k]+X[nodeNum(2)][k]);
    edgeCnt[2][k] = 0.5*(X[nodeNum(0)][k]+X[nodeNum(3)][k]);
    edgeCnt[3][k] = 0.5*(X[nodeNum(1)][k]+X[nodeNum(2)][k]);
    edgeCnt[4][k] = 0.5*(X[nodeNum(1)][k]+X[nodeNum(3)][k]);
    edgeCnt[5][k] = 0.5*(X[nodeNum(2)][k]+X[nodeNum(3)][k]);
  }

  for (int i = 0; i < 4; ++i) {
    int oppi = 3 - i;
    int oppn[3]  = {faceDefTet[oppi][0],faceDefTet[oppi][1],faceDefTet[oppi][2]};
    oppn[0]=nodeNum(oppn[0]); oppn[1]=nodeNum(oppn[1]); oppn[2]=nodeNum(oppn[2]);
    for (int k = 0; k < 3; ++k)
      faceCnt[i][k] = 1.0/3.0*(X[ oppn[0] ][k] + X[ oppn[1] ][k]+X[ oppn[2] ][k]);
  }

  Vec3D hexNodes[8];
  int map1[4] = {0,3,5,2};
  int map2[4] = {3, 3, 3, 1};
  int map3[4] = {1, 0, 3, 5};
  int map4[4] = {2, 4, 5, 4};
  int map5[4] = {2, 0, 1, 0};
  int map6[4] = {1, 2, 0, 2};
  double res[dim];
  
  for (int i = 0; i < 4; ++i) {

    hexNodes[0] = X[ nodeNum(i) ];
    hexNodes[1] = edgeCnt[ map1[i] ];
    hexNodes[2] = faceCnt[ map2[i] ];
    hexNodes[3] = edgeCnt[ map3[i] ];

    
    hexNodes[4] = edgeCnt[ map4[i] ];
    hexNodes[5] = faceCnt[ map5[i] ];
    hexNodes[6] = centroid;
    hexNodes[7] = faceCnt[ map6[i] ];
    
    double eta[3];
    Vec3D xyz;
    for (int j = 0; j < npt; ++j) {
      eta[0] = locs[npt-1][j]*0.5+0.5;
      for (int k = 0; k < npt; ++k) {
	eta[1] = locs[npt-1][k]*0.5+0.5;
	for (int l = 0; l < npt; ++l) {
	  eta[2] = locs[npt-1][l]*0.5+0.5;
	  xyz = hexNodes[0]*(1.0-eta[0])*(1.0-eta[1])*(1.0-eta[2]) +  
	    hexNodes[1]*eta[0]*(1.0-eta[1])*(1.0-eta[2]) +
	    hexNodes[2]*eta[0]*eta[1]*(1.0-eta[2]) +
	    hexNodes[3]*(1.0-eta[0])*eta[1]*(1.0-eta[2]) +

	    hexNodes[4]*(1.0-eta[0])*(1.0-eta[1])*eta[2] +  
	    hexNodes[5]*eta[0]*(1.0-eta[1])*eta[2] +
	    hexNodes[6]*eta[0]*eta[1]*eta[2] +
	    hexNodes[7]*(1.0-eta[0])*eta[1]*eta[2];
	  
	  (obj->*F)( nodeNum(i), xyz, res);
	  for (int m = 0; m < dim; ++m) {
	    V[ nodeNum(i) ][m] += 0.25*vol*wgts[npt-1][j]*wgts[npt-1][k]*wgts[npt-1][l]/8.0*res[m];
	  }
	}
      }
    }
  }
}
