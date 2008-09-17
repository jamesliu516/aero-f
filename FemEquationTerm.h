#ifndef _FEM_EQUATION_TERM_H_
#define _FEM_EQUATION_TERM_H_

#include <BcDef.h>
#include <WallFcn.h>
#include <Vector.h>

// Included
class Communicator;

struct Vec3D;

//------------------------------------------------------------------------------

class FemEquationTerm {

// Included (MB)
public:
  bool completeJac;

protected:

  WallFcn* wallFcn;
  map<int, PorousMedia *> volInfo;

public:

  FemEquationTerm(map<int, VolumeData *> &volData) {
    wallFcn = 0; 
    //construction of volInfo (map from id to porousMedia)
    //...loop on all the VolumeData and check which one is a PorousMedia...
    map<int, VolumeData *>::iterator it;
    if(!volData.empty()){
      for(it=volData.begin(); it!=volData.end(); it++){
        //...if it is a PorousMedia, add it to volInfo
        if(it->second->type == VolumeData::POROUS){
          //...check if it already exists...
          map<int, PorousMedia *>::iterator pmit = volInfo.find(it->first);
          //...if it does not exist yet, add it...
          if(pmit == volInfo.end()) 
            volInfo[it->first] = &(it->second->porousMedia);
        }
      
      }
    }
  }

  ~FemEquationTerm() { if (wallFcn) delete wallFcn; }

  virtual double computeViscousTimeStep(double *, double *) = 0;

  virtual bool computeVolumeTerm(double dp1dxj[4][3], double d2w[4], double *v[4],
				 double *r, double *s, double *, double, SVec<double,3> &, int [4], int) = 0;
  virtual bool computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], double *v[4], 
					 double *drdu, double *dsdu, double *dpdu, double, SVec<double,3> &, int [4], int) = 0;
  virtual void computeSurfaceTerm(int c, Vec3D &n, double d2w[3], 
				  double *vw, double *v[3], double *r) = 0;
  virtual void computeJacobianSurfaceTerm(int c, Vec3D &n, double d2w[3], 
					  double *vw, double *v[3], double *drdu) = 0;
  virtual void computeSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
				  double *vw, double *v[4], double *r) = 0;
  virtual void computeJacobianSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
					  double *vw, double *v[4], double *drdu) = 0;
  virtual bool doesFaceTermExist(int code) {
    if (wallFcn)
      return (code == BC_ADIABATIC_WALL_MOVING || 
	      code == BC_ADIABATIC_WALL_FIXED ||
	      code == BC_ISOTHERMAL_WALL_MOVING || 
	      code == BC_ISOTHERMAL_WALL_FIXED) ? true : false; 
    else
      return (code == BC_ADIABATIC_WALL_MOVING || 
	      code == BC_ADIABATIC_WALL_FIXED) ? true : false; 
  }
  virtual bool doesFaceNeedGradientP1Function() {
    if (wallFcn) return false;
    else return true;
  }
  virtual bool doesSourceTermExist() { return false; }

// Included (MB)
  virtual bool computeDerivativeOfVolumeTerm(double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4], double *v[4],
				 double *dv[4], double dMach, double *dr, double *ds, double *dpr, double dtetvol, SVec<double,3> &x, int nodenum[4], int volid) = 0;
  virtual void computeDerivativeOfSurfaceTerm(int c, Vec3D &n, Vec3D &dn, double d2w[3],
				  double *vw, double *dvw, double *v[3], double *dv[3], double dMach, double *dr) = 0;
  virtual void computeDerivativeOfSurfaceTerm(double dp1dxj[4][3], double ddp1dxj[4][3], int c, Vec3D &n, Vec3D &dn, double d2w[4],
				  double *vw, double *dvw, double *v[4], double *dv[4], double dMach, double *dr) = 0;
  virtual void rstVar(IoData &ioData, Communicator *com) = 0;
  virtual void computeBCsJacobianWallValues(int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]) = 0;

  virtual double computeDerivativeOfViscousTimeStep(double *, double *, double *, double *, double) = 0;


  void computePermittivityTensor(double alpha[3], double beta[3], double ucg[3], double R[9], double *K)
  {
   
    double onehalf = 1.0/2.0;
    double norm_u = pow(ucg[0]*ucg[0] + ucg[1]*ucg[1] + ucg[2]*ucg[2], onehalf); 
    
    double diag[3] = {alpha[0]*norm_u + beta[0], 
                      alpha[1]*norm_u + beta[1],
                      alpha[2]*norm_u + beta[2]};

    double Rt_diag[9] = { R[0]*diag[0], R[3]*diag[1], R[6]*diag[2],
                          R[1]*diag[0], R[4]*diag[1], R[7]*diag[2],
                          R[2]*diag[0], R[5]*diag[1], R[8]*diag[2] };   

    for(int i=0; i<9; ++i)  K[i] = 0.0;
 
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
         for (int k=0; k<3; ++k)
             K[3*i+j] += Rt_diag[3*i+k]*R[3*k+j];

  }

  void computeGradPermittivityTensor(double alpha[3], double ucg[3], double R[9] ,double *B)
  {

    double onehalf = 1.0/2.0;
    double onefourth = 1.0/4.0;

    double norm_u = pow(ucg[0]*ucg[0] + ucg[1]*ucg[1] + ucg[2]*ucg[2], onehalf);
    double inv_norm = 1.0/norm_u;
                                                                                                                                                         
    double diag[3] = {onefourth*inv_norm*alpha[0], 
                      onefourth*inv_norm*alpha[1], 
                      onefourth*inv_norm*alpha[2]};
                                                                                                                                                         
    double Rt_diag[9] = { R[0]*diag[0], R[3]*diag[1], R[6]*diag[2],
                          R[1]*diag[0], R[4]*diag[1], R[7]*diag[2],
                          R[2]*diag[0], R[5]*diag[1], R[8]*diag[2] };
                                                                                                                                                         
    for(int i=0; i<9; ++i)  B[i] = 0.0;
                                                                                                                                                         
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
         for (int k=0; k<3; ++k)
             B[3*i+j] += Rt_diag[3*i+k]*R[3*k+j];

  }


  void multiplyBydVdU(double V[4], double mat[9], double* res, double scaling) 
  {

   double rho = V[0];
   double invRho = 1.0/rho;
   double u = V[1];
   double v = V[2];
   double w = V[3];

   res[0] = 0.0;
   res[1] = 0.0;
   res[2] = 0.0;
   res[3] = 0.0;
   res[4] = 0.0;

   res[5] = -invRho*(mat[0]*u + mat[1]*v + mat[2]*w)*scaling;
   res[6] = mat[0]*invRho*scaling;
   res[7] = mat[1]*invRho*scaling;
   res[8] = mat[2]*invRho*scaling;
   res[9] = 0.0;

   res[10] = -invRho*(mat[3]*u + mat[4]*v + mat[5]*w)*scaling;
   res[11] = mat[3]*invRho*scaling;
   res[12] = mat[4]*invRho*scaling;
   res[13] = mat[5]*invRho*scaling;
   res[14] = 0.0;

   res[15] = -invRho*(mat[6]*u + mat[7]*v + mat[8]*w)*scaling;
   res[16] = mat[6]*invRho*scaling;
   res[17] = mat[7]*invRho*scaling;
   res[18] = mat[8]*invRho*scaling;
   res[19] = 0.0;

   res[20] = 0.0;
   res[21] = 0.0;
   res[22] = 0.0;
   res[23] = 0.0;
   res[24] = 0.0;

  }

/*
  void computeTranformationMatrix(Vec3D iprime, Vec3D jprime, Vec3D kprime, double *R)
  {
     Vec3D nx, ny, nz;

     nx.v[0] = 1.0; nx.v[1] = 0.0; nx.v[2] = 0.0;
     ny.v[0] = 0.0; ny.v[1] = 1.0; ny.v[2] = 0.0;
     nz.v[0] = 0.0; nz.v[1] = 0.0; nz.v[2] = 1.0;

     R[0] = nx*iprime; R[1] = ny*iprime; R[2] = nz*iprime;
     R[3] = nx*jprime; R[4] = ny*jprime; R[5] = nz*jprime;
     R[6] = nx*kprime; R[7] = ny*kprime; R[8] = nz*kprime;

    // to be used in pre-processing stage

     double onehalf = 1.0/2.0;
     Vec3D ny, nz, jprime, kprimei, temp1, temp2;
     ny.v[0] = 0.0; ny.v[1] = 1.0; ny.v[2] = 0.0;
     nz.v[0] = 0.0; nz.v[1] = 0.0; nz.v[2] = 1.0;
     
     temp1 = ny^iprime;
     double v1 = pow(temp1.v[0]*temp1.v[0] + temp1.v[1]*temp1.v[1] + temp3.v[2]*temp3.v[2], onehalf);

     temp2 = nz^iprime;
     double v2 = pow(temp2.v[0]*temp2.v[0] + temp2.v[1]*temp2.v[1] + temp3.v[2]*temp3.v[2], onehalf);
     
     if (v1 > v2) jprime = temp1;
     else jprime = temp2;

     kprime = iprime^jprime;
   
     R[0][0] = 

  }
*/      

};

//------------------------------------------------------------------------------

#endif
