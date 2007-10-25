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
#include <math.h>
#include <GeoState.h>


//------------------------------------------------------------------------------
//--------------functions in ElemTet class
//------------------------------------------------------------------------------

template<int dim>
void ElemTet::computeGalerkinTerm(FemEquationTerm *fet, SVec<double,3> &X, 
				  Vec<double> &d2wall, SVec<double,dim> &V, 
				  SVec<double,dim> &R)
{
  
  double dp1dxj[4][3];
  double vol = computeGradientP1Function(X, dp1dxj);

  double d2w[4] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)],
                   d2wall[nodeNum(2)], d2wall[nodeNum(3)]};
  double *v[4] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)], V[nodeNum(3)]};

  double r[3][dim], s[dim], pr[12];
  bool porousTermExists =  fet->computeVolumeTerm(dp1dxj, d2w, v, reinterpret_cast<double *>(r),
                                                  s, pr, vol, X, nodeNum(), volume_id);

  for (int j=0; j<4; ++j) {
    int idx = nodeNum(j);
    for (int k=0; k<dim; ++k)
      R[idx][k] += vol * ( (r[0][k] * dp1dxj[j][0] + r[1][k] * dp1dxj[j][1] +
                            r[2][k] * dp1dxj[j][2]) - fourth * s[k] );
  }

  if (porousTermExists) {
    for (int j=0; j<4; ++j) {
      int idx = nodeNum(j);
      for (int k=1; k<4; ++k)
        R[idx][k] += pr[3*j+k-1];
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
void ElemTet::computeP1Avg(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test, SVec<double,6> &Eng_Test, 
			   SVec<double,3> &X, SVec<double,dim> &V, double gam, double R)

{

  double dp1dxj[4][3], dudxj[4][3], u[4][3];
  double NCG;
  double Int;                         // stores intermediate values
  double gam1 = gam - 1.0;
  double vol = computeGradientP1Function(X, dp1dxj);

  int i,j,k,l;                        // incrementers for the loops

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

  int i_mom, i_eng;
  i_mom = 0; i_eng = 0;


  // incrementing the counter that stores the volume sum //

  for (i=0; i<4; ++i)
    Mom_Test[nodeNum(i)][i_mom] += vol;

  i_mom =1;

  // adding rho_u computed at the cg to each node //

  for (i=1; i<4; ++i){
    Int = 0.0;
    for (j=0; j<4; ++j){
      Int += V[nodeNum(j)][0]*V[nodeNum(j)][i];
    }
    for (j=0; j<4; ++j){
      Mom_Test[nodeNum(j)][i_mom+i-1] += (NCG*Int);
    }
  }

  i_mom =3;

  // adding rho_u_u computed at the cg to each node //

  l = 0;
  for (i=1; i<4; ++i){
    for (k=i; k<4; ++k){
      Int = 0.0;
      for (j=0; j<4; ++j){
        Int += V[nodeNum(j)][0]*V[nodeNum(j)][i]*V[nodeNum(j)][k];
      }
      l+=1;
      for (j=0; j<4; ++j){
        Mom_Test[nodeNum(j)][i_mom+l] += (NCG*Int);
      }
    }
  }

  i_mom = 10;

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
    for (j=0; j<4; ++j){
      Mom_Test[nodeNum(j)][i_mom+i] += (NCG*Int);
    }
  }

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
  
  i_eng = 1;

  // adding rho_e_plus_p computed at the cg to each node //

  Int = 0.0;
  for (j=0; j<4; ++j){
    Int += V[nodeNum(j)][0]*((1.0/gam1)*(V[nodeNum(j)][4]/V[nodeNum(j)][0])+0.5*squ[j])+V[nodeNum(j)][4];
  }
  for (j=0; j<4; ++j){
    Eng_Test[nodeNum(j)][i_eng] += (NCG*Int);
  }

  i_eng = 2;

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
    for (j=0; j<4; ++j){
      Eng_Test[nodeNum(j)][i_eng+i] += (NCG*Int);
    }
  }

  // incrementing the counter that stores the number of tetrahedrons surrounding a node //

  for (i=0; i<4; ++i)
    Eng_Test[nodeNum(i)][5] += 1.0;

}

//------------------------------------------------------------------------------

template<int dim>
void ElemTet::computeCsValues(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test,
			      SVec<double,6> &Eng_Test, SVec<double,2> &Cs, Vec<double> &VolSum,
			      SVec<double,3> &X, double gam, double R)

{

  double dp1dxj[4][3], dudxj[3][3], u[4][3], ucg[3];
  double vol = computeGradientP1Function(X, dp1dxj);
  double *VC[4] = {VCap[nodeNum(0)], VCap[nodeNum(1)], VCap[nodeNum(2)], VCap[nodeNum(3)]};
  double *mom_test[4] = {Mom_Test[nodeNum(0)], Mom_Test[nodeNum(1)], Mom_Test[nodeNum(2)], Mom_Test[nodeNum(3)]};
  double *eng_test[4] = {Eng_Test[nodeNum(0)], Eng_Test[nodeNum(1)], Eng_Test[nodeNum(2)], Eng_Test[nodeNum(3)]};
  double ntet[4] = {Mom_Test[nodeNum(0)][0], Mom_Test[nodeNum(1)][0], Mom_Test[nodeNum(2)][0], Mom_Test[nodeNum(3)][0]};

// steps in dynamic les to compute unknown coefficients //

  double vc[5];
  double r_u[3];
  double r_u_u[6];
  double r_s_p[6];
  double sq_rat_delta[4] = {pow(eng_test[0][5],(2.0/3.0)), pow(eng_test[1][5],(2.0/3.0)),
                            pow(eng_test[2][5],(2.0/3.0)), pow(eng_test[3][5],(2.0/3.0))};
  double ratdelta;
  double num, denom;
  double sqrt2S2;
  double Pij[3][3], Bij[3][3], Lij[3][3];
  double cs[4], pt[4];

  computeVelocity(VC, u, ucg, ntet);
  computeVelocityGradient(dp1dxj, u, dudxj);
  sqrt2S2 = computeNormSij(dudxj);
  computePij(dudxj, Pij);

// Compute Smagorinsky Coefficient at each node
// ----------------------------------------------
// Ref: Large Eddy Simulation of Bluff-Body flow on Unstructured Grids
// International Journal of Numerical Methods in Fluids 2002
// Vol : 40, pgs:1431-1460
// Authors; Camarri, Salvetti, Koobus, Dervieux
// ---------------------------------------------------------------------


// computing smagorinsky coefficient //

   for (int i=0; i<4; ++i){
     r_u[0] = mom_test[i][1]/mom_test[i][0];
     r_u[1] = mom_test[i][2]/mom_test[i][0];
     r_u[2] = mom_test[i][3]/mom_test[i][0];

     r_u_u[0] = mom_test[i][4]/mom_test[i][0];
     r_u_u[1] = mom_test[i][5]/mom_test[i][0];
     r_u_u[2] = mom_test[i][6]/mom_test[i][0];
     r_u_u[3] = mom_test[i][7]/mom_test[i][0];
     r_u_u[4] = mom_test[i][8]/mom_test[i][0];
     r_u_u[5] = mom_test[i][9]/mom_test[i][0];

     r_s_p[0] = mom_test[i][10]/mom_test[i][0];
     r_s_p[1] = mom_test[i][11]/mom_test[i][0];
     r_s_p[2] = mom_test[i][12]/mom_test[i][0];
     r_s_p[3] = mom_test[i][13]/mom_test[i][0];
     r_s_p[4] = mom_test[i][14]/mom_test[i][0];
     r_s_p[5] = mom_test[i][15]/mom_test[i][0];

     vc[0] = VC[i][0]/mom_test[i][0];
     vc[1] = VC[i][1]/mom_test[i][0];
     vc[2] = VC[i][2]/mom_test[i][0];
     vc[3] = VC[i][3]/mom_test[i][0];
     vc[4] = VC[i][4]/mom_test[i][0];

     ratdelta = sq_rat_delta[i];

     computeLij(Lij, r_u, r_u_u, vc);
     computeBij(Bij, r_s_p, sqrt2S2, Pij, ratdelta, vc);

     num = 0.0;
     denom = 0.0;

     for (int j=0; j<3; ++j){
       for (int k=0; k<3; ++k){
          num += Bij[j][k]*Lij[j][k];
          denom += Bij[j][k]*Bij[j][k];
        }
     }

     if(denom < 0.0000001) denom = 0.0000001;
     cs[i] = (num/denom);

   }


// computing the subgrid scale prandtl number //

   double r_e, r_e_plus_p;
   double r_s_dtdxj[3], dtdxj[3];
   double Li[3], Zi[3];
   double t[4];

   computeTemp(VC, t, R);
   computeTempGradient(dp1dxj, t, dtdxj);


   for (int i=0; i<4; ++i){
     r_u[0] = mom_test[i][1]/mom_test[i][0];
     r_u[1] = mom_test[i][2]/mom_test[i][0];
     r_u[2] = mom_test[i][3]/mom_test[i][0];

     vc[0] = VC[i][0]/mom_test[i][0];
     vc[1] = VC[i][1]/mom_test[i][0];
     vc[2] = VC[i][2]/mom_test[i][0];
     vc[3] = VC[i][3]/mom_test[i][0];
     vc[4] = VC[i][4]/mom_test[i][0];

     r_e = eng_test[i][0]/mom_test[i][0];
     r_e_plus_p = eng_test[i][1]/mom_test[i][0];

     r_s_dtdxj[0] = eng_test[i][2]/mom_test[i][0];
     r_s_dtdxj[1] = eng_test[i][3]/mom_test[i][0];
     r_s_dtdxj[2] = eng_test[i][4]/mom_test[i][0];

     ratdelta = sq_rat_delta[i];

     computeLi(Li, r_e, r_e_plus_p, r_u, vc);
     computeZi(Zi, ratdelta, sqrt2S2, dtdxj, r_s_dtdxj, vc, gam, R);

     num = 0.0;
     denom = 0.0;

     for (int j=0; j<3; ++j){
       num += Li[j]*Zi[j];
       denom += Li[j]*Li[j];
     }

     if(denom < 0.0000001) denom = 0.0000001;
     pt[i] = (num/denom); 

  }

  //----------------------------------------
  // averaging the value of Cs (smoothing)
  //----------------------------------------

  double Int1, Int2;


  Int1 = 0.25*(cs[0]+cs[1]+cs[2]+cs[3]); // avg in each tetrahedron
  Int2 = 0.25*(pt[0]+pt[1]+pt[2]+pt[3]); // avg in each tetrahedron

  for(int i=0; i<4; ++i) {
     Cs[nodeNum(i)][0] += vol*Int1;
     Cs[nodeNum(i)][1] += vol*Int2;
  }

  // incrementing the counter that stores the volume sum

  for (int i=0; i<4; ++i)
    VolSum[nodeNum(i)] += vol;

}

//-----------------------------------------------------------------------

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
void ElemTet::computeDynamicLESTerm(DynamicLESTerm *dles, SVec<double,2> &Cs, Vec<double> &VolSum,
				    SVec<double,3> &X, SVec<double,dim> &V, SVec<double,dim> &R)
{

  double dp1dxj[4][3];
  double vol = computeGradientP1Function(X, dp1dxj);
  double *v[4] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)], V[nodeNum(3)]};
  double cs[4] = {Cs[nodeNum(0)][0]/VolSum[nodeNum(0)], Cs[nodeNum(1)][0]/VolSum[nodeNum(1)],
                  Cs[nodeNum(2)][0]/VolSum[nodeNum(2)], Cs[nodeNum(3)][0]/VolSum[nodeNum(3)]};
  double pt[4] = {Cs[nodeNum(0)][1]/VolSum[nodeNum(0)], Cs[nodeNum(1)][1]/VolSum[nodeNum(1)],
                  Cs[nodeNum(2)][1]/VolSum[nodeNum(2)], Cs[nodeNum(3)][1]/VolSum[nodeNum(3)]};

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
					  SVec<double,dim> &V, GenMat<Scalar,neq> &A)
{

  double dp1dxj[4][3];
  double vol = computeGradientP1Function(X, dp1dxj);

  double d2w[4] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)],
                   d2wall[nodeNum(2)], d2wall[nodeNum(3)]};
  double *v[4] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)], V[nodeNum(3)]};

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
