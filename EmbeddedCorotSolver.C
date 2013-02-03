#include <cstdlib>
#include <cmath>

#include <EmbeddedCorotSolver.h>

#include <MatchNode.h>
#include <Domain.h>
#include <Vector3D.h>
#include <DenseMatrixOps.h>

#include <BCApplier.h>

//------------------------------------------------------------------------------

EmbeddedCorotSolver::EmbeddedCorotSolver(DefoMeshMotionData &data, Domain *dom, double *Xstruct, int nNodes)
  : domain(dom), Xs0(Xstruct), X0(dom->getNodeDistInfo())
{

  numStNodes = nNodes;

  numLocSub = domain->getNumLocSub();

  com = domain->getCommunicator();

  domain->getReferenceMeshPosition(X0);

  computeCG(Xs0, cg0);

  double zeroRot[3] = {0.0, 0.0, 0.0};
  computeRotMat(zeroRot, R);

  //HB: look if a symmetry plane was specified in the input file
  double nx = data.symmetry.nx;
  double ny = data.symmetry.ny;
  double nz = data.symmetry.nz;
  double nrm= sqrt(nx*nx+ny*ny+nz*nz);
  SymAxis = EmbeddedCorotSolver::NONE;
  if(nrm!=0.0){
    nx /= nrm; ny /= nrm; nz /= nrm; 
    if((fabs(nx)==1.0) & (ny==0.0) & (nz==0.0))
      SymAxis = EmbeddedCorotSolver::AXIS_X;
    else if((nx==0.0) & (fabs(ny)==1.0) & (nz==0.0))
      SymAxis = EmbeddedCorotSolver::AXIS_Y;
    else if((nx==0.0) & (ny==0.0) & (fabs(nz)==1.0))
      SymAxis = EmbeddedCorotSolver::AXIS_Z;
    else {
      com->fprintf(stderr," *** ERROR: embedded corotational solver only supports a canonical plane as a symmetry plane.\n");
      exit(-1);
    }
  } 
  switch(SymAxis) {
    case(EmbeddedCorotSolver::NONE):
      com->fprintf(stderr," ... No symmetry plane is used in the embedded corotational solver.\n");
      com->fprintf(stderr,"     -> the 3 rotations axis (X,Y,Z) & 3 translation axis (X,Y,Z) are used.\n");
      break;
    case(EmbeddedCorotSolver::AXIS_X):
      com->fprintf(stderr," ... Symmetry plane of normal X is used in the embedded corotational solver.\n");
      com->fprintf(stderr,"     -> only rotation around axis X & translations in the Y-Z plane are allowed.\n");
      break;
    case(EmbeddedCorotSolver::AXIS_Y):
      com->fprintf(stderr," ... Symmetry plane of normal Y is used in the embedded corotational solver.\n");
      com->fprintf(stderr,"     -> only rotation around axis Y & translations in the X-Z plane are allowed.\n");
      break;
    case(EmbeddedCorotSolver::AXIS_Z):
      com->fprintf(stderr," ... Symmetry plane of normal Z is used in the embedded corotational solver.\n");
      com->fprintf(stderr,"     -> only rotation around axis Z & translations in the X-Y plane are allowed.\n");
      break;
  }
}

//------------------------------------------------------------------------------

inline
void invRotLocVec(double mat[3][3], double v[3]) 
{

  double c[3];
  for (int j = 0; j < 3; j++)
    c[j] = mat[0][j]*v[0] + mat[1][j]*v[1] + mat[2][j]*v[2];

  v[0] = c[0];
  v[1] = c[1];
  v[2] = c[2];

}

//------------------------------------------------------------------------------

void EmbeddedCorotSolver::computeRotGradAndJac(double *Xs, double RR[3][3], 
					 double cg1[3], double grad[3], double jac[3][3])
{
  int i, j;

  for (i = 0; i < 3; i++)  {
    grad[i] = 0.0;
    for (j = 0; j < 3; j++)
      jac[i][j] = 0.0;
  }

  Vec3D rotGrad[3];

  for (i = 0; i < numStNodes; i++) {
    double rd[3];
    // rotate the local vectors using R(n-1)
    rd[0] = Xs0[3*i+0] - cg0[0];
    rd[1] = Xs0[3*i+1] - cg0[1];
    rd[2] = Xs0[3*i+2] - cg0[2];
    rotLocVec(RR, rd);

    Vec3D eVec;
    //compute freq. used values
    eVec[0] = Xs[3*i+0] - cg1[0] - rd[0];
    eVec[1] = Xs[3*i+1] - cg1[1] - rd[1];
    eVec[2] = Xs[3*i+2] - cg1[2] - rd[2];

    rotGrad[0][0] = -rd[1];
    rotGrad[0][1] = rd[0];
    rotGrad[0][2] = 0;
  
    rotGrad[1][0] = rd[2];
    rotGrad[1][1] = 0;
    rotGrad[1][2] = -rd[0];

    rotGrad[2][0] = 0;
    rotGrad[2][1] = -rd[2];
    rotGrad[2][2] = rd[1];

    grad[0] += -2*(eVec * rotGrad[0]);
    grad[1] += -2*(eVec * rotGrad[1]);
    grad[2] += -2*(eVec * rotGrad[2]);

    jac[0][0] +=  2 * (rd[0]*eVec[0]
                    +  rd[1]*eVec[1]
                    +  rotGrad[0][0]*rotGrad[0][0]
                    +  rotGrad[0][1]*rotGrad[0][1]);

    jac[0][1] += -2 * (rd[2]*eVec[1]
                    - rotGrad[0][0]*rotGrad[1][0]);

    jac[0][2] += -2 * (rd[2]*eVec[0]
                    -  rotGrad[0][1]*rotGrad[2][1]);

    jac[1][1] += 2 * (rd[0]*eVec[0]
                    +  rd[2]*eVec[2]
                    +  rotGrad[1][0]*rotGrad[1][0]
                    +  rotGrad[1][2]*rotGrad[1][2]);

    jac[1][2] += -2 * (rd[1]*eVec[0]
                    -   rotGrad[2][2]*rotGrad[1][2]);

    jac[2][2] += 2 * (rd[2]*eVec[2]
                    +  rd[1]*eVec[1]
                    +  rotGrad[2][1]*rotGrad[2][1]
                    +  rotGrad[2][2]*rotGrad[2][2]);
  }

  // Fill in symmetric terms
  jac[1][0] = jac[0][1];
  jac[2][0] = jac[0][2];
  jac[2][1] = jac[1][2];

}

//------------------------------------------------------------------------------

void EmbeddedCorotSolver::rotLocVec(double mat[3][3], double v[3]) 
{
  double c[3];
  for (int j = 0; j < 3; j++)
    c[j] = mat[j][0]*v[0] +
           mat[j][1]*v[1] +
           mat[j][2]*v[2];
  
  v[0] = c[0];
  v[1] = c[1];
  v[2] = c[2];
}

//------------------------------------------------------------------------------
//computes delta R

void EmbeddedCorotSolver::computeRotMat(double *angle, double mat[3][3])
{
  // trig functions of angles
  double c1 = cos(angle[0]);
  double s1 = sin(angle[0]);
  double c2 = cos(angle[1]);
  double s2 = sin(angle[1]);
  double c3 = cos(angle[2]);
  double s3 = sin(angle[2]);

  /* 
     compute rotation matrix computed as R1.R2.R3
     where R1 is rotation about z
           R2 is rotation about y
           R3 is rotation about x
  */

  mat[0][0] = c1*c2;
  mat[0][1] = c1*s2*s3 - c3*s1;
  mat[0][2] = c1*c3*s2 + s1*s3;

  mat[1][0] = c2*s1;
  mat[1][1] = c1*c3+s1*s2*s3;
  mat[1][2] = c3*s1*s2 - c1*s3;

  mat[2][0] = -s2;
  mat[2][1] = c2*s3;
  mat[2][2] = c2*c3;
}

//------------------------------------------------------------------------------

void EmbeddedCorotSolver::printRotMat(double mat[3][3])
{
  com->fprintf(stderr," Rotation matrix = \n");
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++)
      com->fprintf(stderr," %e  ",mat[i][j]);
    com->fprintf(stderr,"\n");
  }
}

//------------------------------------------------------------------------------

void EmbeddedCorotSolver::computeCG(double *Xs, double cg[3])
{
  
  for (int j=0; j<3; ++j) cg[j] = 0.0;

  for (int i=0; i<numStNodes; ++i) {
    for (int j=0; j<3; ++j)
      cg[j] += Xs[3*i+j];
  }

  double invTotNd = 1.0 / double(numStNodes);

  for (int j=0; j<3; ++j)
    cg[j] *= invTotNd;
}

//------------------------------------------------------------------------------

void EmbeddedCorotSolver::computeNodeRot(double RR[3][3], DistSVec<double,3> &X, 
				 double cg00[3], double cg1[3])
{
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)  {

    double (*x)[3]  = X.subData(iSub);
    double (*x0)[3] = X0.subData(iSub);

    for (int i = 0; i < X.subSize(iSub); ++i)  {
      x[i][0] = x0[i][0] - cg00[0];
      x[i][1] = x0[i][1] - cg00[1];
      x[i][2] = x0[i][2] - cg00[2];

      rotLocVec(RR, x[i]);

      x[i][0] = x[i][0] + cg1[0];
      x[i][1] = x[i][1] + cg1[1];
      x[i][2] = x[i][2] + cg1[2];
    }
  }
}

//------------------------------------------------------------------------------

void EmbeddedCorotSolver::solveDeltaRot(double *Xs, double cg1[3])
{
  double jac[3][3], grad[3];
  double dRot[3][3];

  int maxits = 10;
  double atol = 1.e-12;
  double rtol = 1.e-10;
  double res0, res, target;

  for (int iter = 0; iter < maxits; ++iter) {

    // compute rotation gradients and derivatives
    computeRotGradAndJac(Xs, R, cg1, grad, jac);

    //HB: zero terms depending on the axis of the plane of symmetry
    switch(SymAxis) {
      case(EmbeddedCorotSolver::AXIS_X): 
        grad[1] = grad[2] = 0.0;
        jac[0][1] = jac[1][0] = jac[0][2] = jac[2][0] = 0.0;
        jac[1][2] = jac[2][1] = 0.0;
        jac[1][1] = jac[2][2] = 1.0;
        break;
      
      case(EmbeddedCorotSolver::AXIS_Y):
        grad[0] = grad[2] = 0.0;
        jac[0][1] = jac[1][0] = jac[0][2] = jac[2][0] = 0.0;
        jac[1][2] = jac[2][1] = 0.0;
        jac[0][0] = jac[2][2] = 1.0;
        break;

      case(EmbeddedCorotSolver::AXIS_Z):
        grad[0] = grad[1] = 0.0;
        jac[0][1] = jac[1][0] = jac[0][2] = jac[2][0] = 0.0;
        jac[1][2] = jac[2][1] = 0.0;
        jac[0][0] = jac[1][1] = 1.0;
        break;
    }

    res = sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);

    if (iter == 0) { res0 = res; target = rtol*res0; } 
    
    //if (res == 0.0 || res <= target) break;
    if (res <= atol || res <= target) break; //HB: add absolute tol. to avoid iterations
                                             //if the initial residual is already small

    // rotation results come back in grad
    solveRotMat(jac, grad);

    grad[0] *= -1.0;
    grad[1] *= -1.0;
    grad[2] *= -1.0;

    // update dRot
    computeRotMat(grad, dRot);

    denseMatrixTimesDenseMatrixInPlace(dRot, R);

  }

  if (res>target & res>atol) {
    com->printf(1, "*** Warning: incremental rotation solver reached %d its", maxits);
    com->printf(1, " (initial res = %.2e, final res=%.2e, target=%.2e)\n", res0, res, target);
  }
}

//------------------------------------------------------------------------------

void EmbeddedCorotSolver::solveRotMat(double m[3][3], double v[3])  
{
  int i,j,k;
  for (i = 0; i < 2; i++)
    for (j = i+1; j < 3; ++j) {
      double coef = m[j][i]/m[i][i];
      for (k = i+1; k < 3; ++k)
        m[j][k] -= coef*m[i][k];
      v[j] -= coef*v[i];
    }
  
  for (i=2; i >= 0; i--) {
    for (j=2; j > i; j--)
      v[i] -= m[i][j] * v[j];
    v[i] /= m[i][i];

  }
}

//------------------------------------------------------------------------------

void EmbeddedCorotSolver::solve(double *Xtilde, int nNodes, DistSVec<double,3> &X)
{

  if(nNodes!=numStNodes) {
    com->fprintf(stderr,"Number of structure nodes has changed!\n");
    exit(-1);
  }

  // compute cg(n+1)
  double cg1[3];
  computeCG(Xtilde, cg1);

  switch(SymAxis) {
    case(EmbeddedCorotSolver::AXIS_X):
      cg1[0] = cg0[0];
      break;

    case(EmbeddedCorotSolver::AXIS_Y):
      cg1[1] = cg0[1];
      break;

    case(EmbeddedCorotSolver::AXIS_Z):
      cg1[2] = cg0[2];
      break;
  }

  // solve for the incremental rotations via Newton-Rhapson
  solveDeltaRot(Xtilde, cg1);

// Update node
  computeNodeRot(R, X, cg0, cg1);

}

//------------------------------------------------------------------------------
