/* MultiGridKernel.h

 */

#pragma once

#include <DistMatrix.h>
#include <DistGeoState.h>
#include <MultiGridLevel.h>
#include <KspSolver.h>
#include <DistEmbeddedVector.h>
#include <MatVecProd.h>
#include <MultiGridOperator.h>
#include <DistTimeState.h>
#include <SpaceOperator.h>

class DistGeoState;
template<class Scalar,int dim> class DistSVec;

template<class Scalar,int dim> class MultiGridPrecMatVecProd;
template<class Scalar,int dim,class Scalar2 = double> class MultiGridPrecJacobiPrec;
template<class Scalar,int dim,class Scalar2 = double> class MultiGridPrecRASPrec;

template<class Scalar, int dim,class Scalar2 = double>
class MultiGridKernel {

 public:

  class MultiGridSmoother {

   public:

    MultiGridSmoother(MultiGridKernel<Scalar,dim,Scalar2>* owner) : pOwner(owner) { }

    ~MultiGridSmoother() { }

    virtual void smooth(int level, DistSVec<Scalar2,dim>& x,
                        const DistSVec<Scalar2,dim>& f,int steps) {

      pOwner->smooth(level,x,f,steps);
    }

    virtual void applyOperator(int level, DistSVec<Scalar2,dim>& f,
                             DistSVec<Scalar2,dim>& x) {

      pOwner->applyOperator(level,f,x);
    }

   protected:

    MultiGridKernel<Scalar,dim,Scalar2>* pOwner;
  };

  MultiGridKernel(Domain *dom, DistGeoState& distGeoState, KspData&
                  coarseSolverData,IoData&,VarFcn* varFcn,
                  bool createFineA, int num_levels,
                  DistTimeState<dim>*,BcFcn* bcFcn = NULL);

  ~MultiGridKernel();

  void setParameters(int v1, int v2, int
                     fine_sweeps, double relax, int do_out);

  void setOperators(SpaceOperator<dim>*);

  void initialize();

  void setupAlgebraic();

  void setupGeometric(DistSVec<Scalar2,dim>& U);

  void smooth(int level,DistSVec<Scalar2,dim>& x,
              const DistSVec<Scalar2,dim>& f,int steps);

  void cycleV(DistSVec<Scalar2,dim>& f, 
              DistSVec<Scalar2,dim>& x);
  void cycleW(int lvl,DistSVec<Scalar2,dim>& f, 
              DistSVec<Scalar2,dim>& x);

  void getData(DistMat<Scalar2,dim>& mat);

  DistMat<Scalar2,dim>& getFineMatrix() { return *macroA[0]; }

  bool isInitialized() const { return initialized; }

  void setGeometric();
    
  void applyOperator(int level, DistSVec<Scalar2,dim>& f,
                   DistSVec<Scalar2,dim>& x);

  void setupBcs(DistBcData<dim>*);

  MultiGridOperator<Scalar2,dim>* getOperator(int level) {

    return myOperators[level];
  }

  DistSVec<Scalar2,dim>& getTemporary(int lvl) {

    return *macroValuesTmp[lvl];
  }
  
  DistSVec<Scalar2,dim>& getTemporaryV(int lvl) {

    return *macroValuesTmpV[lvl];
  }
  
  DistSVec<Scalar2,dim>& getResidual(int lvl) {

    return *macroR[lvl];
  }
  
  DistSVec<Scalar2,dim>& getDX(int lvl) {

    return *macroDX[lvl];
  }

  DistMvpMatrix<Scalar2,dim>& getA(int lvl) {

    return dynamic_cast<DistMvpMatrix<Scalar2,dim>&>(*macroA[lvl]);
  }
 
  void kspSolve(int lvl, DistSVec<Scalar2, dim>& f, DistSVec<Scalar2, dim>& x); 

  void setSmoother(MultiGridSmoother* s) { mySmoother = s; }

  void setupPreconditioner(int lvl);

  MultiGridLevel<Scalar2>* getLevel(int lvl) { return multiGridLevels[lvl]; }

  void setProlongRelaxFactor(double r) { prolong_relax_factor = r; }
  void setRestrictRelaxFactor(double r) { restrict_relax_factor = r; }

 private:

  bool isGeometric;

  int nSmooth1,nSmooth2;

  double relaxationFactor;
  
  const int num_levels, agglom_size, numLocSub;

  double beta; 
 
  double prolong_relax_factor;
  double restrict_relax_factor;

  MultiGridLevel<Scalar2> ** multiGridLevels;
 
  MultiGridOperator<Scalar2,dim>** myOperators;
 
  MultiGridSmoothingMatrix<Scalar2,dim> ***smoothingMatrices;

  DistSVec<Scalar2, dim> ** macroValues, ** macroValuesTmp,**macroValuesOld;
  DistSVec<Scalar2, dim> ** macroR;
  DistSVec<Scalar2, dim> ** macroValuesTmpV;
  DistVec<Scalar2> ** macroIrey;
  DistSVec<Scalar2, dim> ** macroDX;

  DistMat<Scalar2,dim> ** macroA;

  DistTimeState<dim>* myTimeState;
  VarFcn* myVarFcn;
  FluxFcn** myFluxFcn;

  DistGeoState& geoState;

  IoData& ioData;

  KspSolver<DistSVec<Scalar2,dim>, 
            MultiGridPrecMatVecProd<Scalar2,dim> ,
            MultiGridPrecRASPrec<Scalar2,dim>, Communicator>** coarseSolvers;

  MultiGridPrecMatVecProd<Scalar2,dim>** coarseMvps;
  MultiGridPrecRASPrec<Scalar2,dim>** coarsePrecs;

  bool ownsFineA;

  Domain* domain;

  int output;

  int fine_sweeps;

  bool initialized;
 
  KspData& coarseSolverData;
  
  MultiGridSmoother* mySmoother;

  MultiGridSmoother defaultSmoother;

  BcFcn* bcFcn;
};
