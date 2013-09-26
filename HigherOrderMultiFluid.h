/* HigherOrderMultiFluid.h

 */

#pragma once

#include <NodalGrad.h>

struct V6NodeData;

class HigherOrderMultiFluid {

  public:

    HigherOrderMultiFluid();

    ~HigherOrderMultiFluid();

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
     double estimateR(int l, int vertex, 
		      int i, SVec<double,dim>& V, 
		      NodalGrad<dim>& dVdx, SVec<double,3>& X,
		      Vec<int>& fluidId);

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
		   Vec<int>& fluidId) ;

  private:

   void* lastPhaseChangeState;

   ElemSet* elems;

   V6NodeData (*v6data)[2];

   bool limitExtrap;
};

#include <HigherOrderMultiFluid.C>
