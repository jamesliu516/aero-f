/* HigherOrderMultiFluid.h

 */

#pragma once

#include <NodalGrad.h>

struct V6NodeData;

class HigherOrderMultiFluid {

  public:

   struct CutCellState {

     int fid1,fid2;
     void* cutCellData;
     
   };

    HigherOrderMultiFluid(Vec<CutCellState*>& myVec);

    ~HigherOrderMultiFluid();
    /*
   void computeLevelSetIntersection(const double x0[3],
                                    const double x1[3],
                                    int subId,
                                    int nodeId,
                                    double iloc[3]); 
  
   void computeExtrapolation(const double x1[3],
                             const double iloc[3],
                             std::pair<int,int> node0,
                             std::pair<int,int> node1,
                             Scalar Uext1[dim],
                             Scalar Uext2[dim]);    

   void computeInterpolation(const double x0[3],
                     const double iloc[3],
                     const double xmid[3],
                     const Scalar U0[dim],
                     const Scalar Ustar[dim],
                     Scalar Ui[dim]);
    */

   int isCellCut(int i) const { return cutCells[i] != 0; }

   int getOtherFluidId(int cell, int fid) {

     CutCellState* C = cutCells[cell];
     return (C->fid1 == fid ? C->fid2 : C->fid1);
   }

   template <int dim>
   void computeCutCellExtrapolations(int cutCellId,int fidi,int fidj, 
				     const double iloc[3],
				     double* Vi, double* Vj,
				     SVec<double,3>& X);

   template <int dim>
     void setCutCellFlags(int lsdim, Vec<int>& status);

   template <int dim>
     void clearCutCellFlags();

   template <int dim>
     void printCutCellData(int i);

   int getNumCutCells();
   
   template <int dim>
     void storeCutCellData(SVec<double,dim>* cutCell[2],
			   NodalGrad<dim,double>* cutGrad[2],
			   Vec<int>* counts[2]);

   template<int dim>
     void setCutCellData(SVec<double,dim>& V, Vec<int>& fid);

   template<int dim>
     void getCutCellData(int,int fid,double V[dim], double x[dim][3]);

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

  private:
   /*
    DistVec<int>* nodeStatus[dimLS];

    DistVec<int> myFid;

    DistSVec<Scalar,dimLS>* levelSet;
    DistSVec<Scalar,dim>* U;

    DistNodalGrad<Scalar,dimLS>* levelSetGrad;

    DistNodalGrad<Scalar,dim>* dU;
   */

   template <int dim>
    struct CutCellStateData {

      // A vector of the extrapolated states at cells that
      // are cut.
      double V[2][dim];
      
      // A vector of gradients at cut cells, to enable extrapolation
      double dV[2][dim][3]; 

    };


   Vec<CutCellState*>& cutCells;

   int numCutCells;

   void* lastPhaseChangeState;

   ElemSet* elems;

   V6NodeData (*v6data)[2];
};

#include <HigherOrderMultiFluid.C>
