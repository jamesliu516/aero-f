/* HigherOrderFSI.h

 */

#pragma once

#include <NodalGrad.h>

struct V6NodeData;

class HigherOrderFSI {

  public:

    HigherOrderFSI();

    ~HigherOrderFSI();

   void setLimitedExtrapolation();

   template<int dim>
     void initialize(int numNodes,ElemSet&, V6NodeData (*)[2]);

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

   //V6NodeData (*getV6Data()) [2] { return v6data; }

   //private:

   void* lastPhaseChangeState;

   ElemSet* elems;

   V6NodeData (*v6data)[2];

   bool limitExtrap;
};

#include <HigherOrderFSI.C>
