#ifndef _HIGHER_ORDER_FSI_H_
#define _HIGHER_ORDER_FSI_H_

#include <NodalGrad.h>
#include <FemEquationTermDesc.h>

struct V6NodeData;

class HigherOrderFSI {

  public:
	int locationhalfriemannproblem;// {INTERSECTPOINT = 0, CLOSESTPOINT = 1} ;


    HigherOrderFSI(const IoData &iod);

    ~HigherOrderFSI();

   void setLimitedExtrapolation();

   template<int dim>
     void initialize(int numNodes,ElemSet&, V6NodeData (*)[2]);

   template <int dim>
		void initialize(IoData &iod, Communicator *comm, int numNodes, ElemSet&, bool linearReconstruction, bool hoViscReconstruction);

   template <int dim>
     bool hasLastPhaseChangeValue(int nodeId);

   template <int dim>
     const double* getLastPhaseChangeValue(int nodeId);

   template <int dim>
     void setLastPhaseChangeValue(int nodeId,const double*);

   template <int dim>
     void setLastPhaseChangeValues(SVec<double,dim>& update,
				   Vec<double>& weight);

   template <int dim>
     void estimateR(int l, int vertex, 
		    int i, SVec<double,dim>& V, 
		    NodalGrad<dim>& dVdx, SVec<double,3>& X,
		    Vec<int>& fluidId, double* r);

   template <int dim>
     void estimateRderivative(int l, int vertex, 
			      int i, SVec<double,dim>& V, 
			      NodalGrad<dim>& dVdx, SVec<double,3>& X,
			      Vec<int>& fluidId, double* dV, double* dV_g);

   template <int dim>
     double computeAlpha(int nodeId, const double* currentV,
			 const double* neighborV);

   bool limitExtrapolation() const { return limitExtrap; }

   template <int dim>
     void
     extrapolateV6(int l, int vertex, 
		   int i, SVec<double,dim>& V, 
		   double* Vsurrogate,const double* W, SVec<double,3>& X,
		   double alpha,double length,
		   Vec<int>& fluidId, double* beta) ;

   template <int dim>
     void RcnExtrap(int l, int vertex, int i, 
		    double length, double alphaij, 
		    SVec<double,dim>& p,
		    double* ddpij,
		    SVec<double,3>& X, Vec<int>& fluidId,
		    double* beta,
		    double* dVsdV, 
		    double *dVpij, double *dVpji);

   
   template <int dim>
     void derivativeofHOFSI(int l, int vertex, int i, 
			    SVec<double,dim>& V, 
			    double*  Vf,    double*  Vs,
			    double* dVf_ds, double* dVs_ds,
			    SVec<double,3>& X,
			    double alphaij, double dalphaij_ds,
			    double length, Vec<int>& fluidId, double* beta,
			    double* dVfluid_ds, double* dVstar_ds);

   V6NodeData (*getV6Data() const) [2] { return v6data; }

	/* --------------------------- */
//
//		void setSIstencil(V6NodeData *SIstencilData);
//
//
//		void setFEMstencil(V6NodeData *FEMstencilData_p,
//								 V6NodeData *FEMstencilData_m);

	template<int dim>
		void extrapolateToWall_1(int l, int n, int Fid, VarFcn *varFun, 
										 SVec<double,dim>& V, NodalGrad<dim>& dV, double* V_n, 
										 SVec<double,3>& X, Vec3D &xWall, Vec3D &nWall, Vec3D &Xij,
										 double* V_ext, bool externalSI);

	template<int dim>
		void interpolateToSI(int l, int n, int Fid, VarFcn *varFun, 
									SVec<double,dim>& V, double* Vstar, NodalGrad<dim>& dV,
									SVec<double,3>& X, Vec3D &xWall, Vec3D &nWall, Vec3D &Xij,
									double* Vsi, double limiter = 0.5, bool externalSI = true);
	template<int dim>
		bool setFEGhostPoint(int dir, int i, VarFcn *varFun, SVec<double,dim>& U, 
									 NodalGrad<dim>& dV, SVec<double,3>& X, Vec<int> &fluidId,
									 Vec3D &xWall, Vec3D &vWall, Vec3D &nWall, 
									 bool isIsoTherm, double TWall, 
									FemEquationTerm *fet, double* Vg, int &fId);

	double vanAlbada(double a, double b);


	void safeExtrapolation(int dim, const double* Vi, const double* Vghost, const double *dV, bool ij, double alpha, double* Ve);

	void safeExtrapolation(int dim, const double* V, const double* dV, double beta, double* Ve);
	/* --------------------------- */

private:

    void* lastPhaseChangeState;

    ElemSet* elems;

    V6NodeData (*v6data)[2];

    V6NodeData (*SIData_p);
	V6NodeData (*SIData_m);
    V6NodeData (*FEMData_p);
    V6NodeData (*FEMData_m);


	double geomTol;
	bool limitExtrap;

	bool HOtreatment;
	bool viscQuadRcn;
public:
	//todo here we assume the node or edge numbers in each subdomain is fixed, need more for AMR
	V6NodeData* getAllocatedSIData_p(int len)  {if(!SIData_p) SIData_p = new V6NodeData[len];  return SIData_p;}
	V6NodeData* getAllocatedSIData_m(int len)  {if(!SIData_m) SIData_m = new V6NodeData[len];  return SIData_m;}
	V6NodeData* getAllocatedFEMData_p(int len) {if(!FEMData_p) FEMData_p = new V6NodeData[len];  return FEMData_p;}
	V6NodeData* getAllocatedFEMData_m(int len) {if(!FEMData_m) FEMData_m = new V6NodeData[len];  return FEMData_m;}


	void setVelocityG(double* Vf, double** dVf,
							Vec3D Dir, double xi, double eta,
							Vec3D &vWall, double* Vg);

	
	void setVelocityWfG(double* Vf, double** dVf, VarFcn *vf,
							  Vec3D Dir, double xi, double eta,
							  Vec3D &vWall, Vec3D &nW, 
							  Vec3D &tgW1, Vec3D &tgW2,
							  double dudn, double* Vg);

	void setTemperatureG(double* Vf, double** dVf,
								Vec3D Dir, double xi, double eta,
								bool isIsoTherm, double TWall, double* Vg);

	void setTurboG(double* Vf, double** dVf,
						Vec3D Dir, double xi, double eta,
						int dim, double* Vg);

	void computeWallVersors(double *V1, Vec3D &nW, VarFcn *vf,
									Vec3D &tgW1, Vec3D &tgW2);

	bool computeDuDTwf(double *V1, VarFcn *vf, double d2w, 
							 Vec3D &vWall, double TWall, 
							 Vec3D &tgW1, Vec3D &tgW2,
							 FemEquationTerm *fet, 
							 double &dudn, double &dTdn);


};

#ifdef TEMPLATE_FIX
#include <HigherOrderFSI.C>
#endif

#endif
