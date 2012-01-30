/* HigherOrderMultiFluid.h

 */

#pragma once

class HigherOrderMultiFluid {

  public:

    HigherOrderMultiFluid();

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

   int getNumCutCells();

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

   struct CutCellState {

     int fid1,fid2;
     void* cutCellData;
     
   };


   Vec<CutCellState*> cutCells;

   int numCutCells;
};

#include <HigherOrderMultiFluid.C>
