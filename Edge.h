#ifndef _EDGE_H_
#define _EDGE_H_

#ifdef OLD_STL
#include <map.h>
#else
#include <map>
using std::map;
#endif

#if defined(__KCC)
#include <utility.h>
#else
#include <utility>
using std::pair;
#endif

// Included (MB)
#include <cstdio>

#include <ProgrammedBurn.h>

class FluidSelector;
class VarFcn;
class RecFcn;
class FluxFcn;
class ElemSet;
class GeoState;
class FemEquationTerm;
class TimeLowMachPrec;
class LevelSetStructure;


struct Vec3D;

#ifndef _NDGRAD_TMPL_
#define _NDGRAD_TMPL_
template<int dim, class Scalar = double> class NodalGrad;
#endif

template<int dim> class EdgeGrad;
template<int dim> class ExactRiemannSolver;
template<class Scalar> class Vec;
template<class Scalar, int dim> class SVec;
template<class Scalar, int dim> class GenMat;

//------------------------------------------------------------------------------

class EdgeSet {

  typedef pair<int, int> Pair;

#ifdef OLD_STL
  typedef map<Pair, int, less<Pair> > MapPair;
#else
  typedef map<Pair, int> MapPair;
#endif

  MapPair *mp;

  int numEdges;

  int (*ptr)[2];
  bool *masterFlag;

  double* edgeLength;

  ProgrammedBurn* programmedBurn;

public:

  EdgeSet();
  ~EdgeSet();

  int find(int, int);

  void createPointers(Vec<int> &);

  template<int dim>
  void computeTimeStep(FemEquationTerm *, VarFcn *, GeoState &,
		       SVec<double,3> &, SVec<double,dim> &, Vec<double> &,
                       Vec<double> &, TimeLowMachPrec &);
  template<int dim>
  void computeTimeStep(VarFcn *, GeoState &, SVec<double,dim> &, Vec<double> &,
                       TimeLowMachPrec &, Vec<int> &, int,Vec<double>*);

  template<int dim>
  int computeFiniteVolumeTerm(int*, Vec<double> &, FluxFcn**, RecFcn*, ElemSet&, GeoState&,
                              SVec<double,3>&, SVec<double,dim>&, NodalGrad<dim>&, EdgeGrad<dim>*,
			      SVec<double,dim>&, SVec<int,2>&, int, int);

  template<int dim, int dimLS>
  int computeFiniteVolumeTerm(ExactRiemannSolver<dim>&, int*,
                              FluxFcn**, RecFcn*, ElemSet&, GeoState&, SVec<double,3>&,
                              SVec<double,dim>&, Vec<int> &, FluidSelector &,
                              NodalGrad<dim>&, EdgeGrad<dim>*,
                              NodalGrad<dimLS>&,
                              SVec<double,dim>&, int,
                              SVec<int,2>&, int, int);

  template<int dim, int dimLS>
  int computeFiniteVolumeTerm(ExactRiemannSolver<dim>&, int*,
                              FluxFcn**, RecFcn*, ElemSet&, GeoState&, SVec<double,3>&,
                              SVec<double,dim>&, SVec<double,dim>&, SVec<double,dim>&,
                              LevelSetStructure&, bool, Vec<int> &, int, SVec<double,3>*, FluidSelector &,
                              NodalGrad<dim>&, EdgeGrad<dim>*,
                              NodalGrad<dimLS>&,
                              SVec<double,dim>&, int,
                              SVec<int,2>&, int, int);

  /** compute flux for Riemann based FSI*/
  template<int dim>
  int computeFiniteVolumeTerm(ExactRiemannSolver<dim>&, int*,
                              FluxFcn**, RecFcn*, ElemSet&, GeoState&, SVec<double,3>&,
                              SVec<double,dim>&, SVec<double,dim>&, SVec<double,dim>&, LevelSetStructure &,
                              bool, Vec<int>&, int, SVec<double,3>*, NodalGrad<dim>&, EdgeGrad<dim>*,
                              SVec<double,dim>&, int,
                              SVec<int,2>&, int, int);

  template<int dim, int dimLS>
  void computeFiniteVolumeTermLS(FluxFcn**, RecFcn*, RecFcn*, ElemSet&, GeoState&, SVec<double,3>&,
                               SVec<double,dim>&, NodalGrad<dim>&, NodalGrad<dimLS>&, EdgeGrad<dim>*,
                               SVec<double,dimLS>&, SVec<double,dimLS>&, LevelSetStructure* =0);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, GeoState &,
                                       Vec<double> &, SVec<double,3> &, Vec<double> &,
                                       SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, GeoState &,
                               Vec<double> &, SVec<double,3> &, Vec<double> &,
                               SVec<double,dim> &, GenMat<Scalar,neq> &,
                               int * );

  template<int dim, class Scalar, int neq, int dimLS>
  void computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>&, FluxFcn **, GeoState &,
                                       NodalGrad<dim> &, NodalGrad<dimLS> &,
                                       SVec<double,3> &, Vec<double> &,
                                       SVec<double,dim> &, GenMat<Scalar,neq> &,
                                       FluidSelector &, Vec<int> &);
  template<int dim, class Scalar, int neq, int dimLS>
  void computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>&, FluxFcn **, GeoState &,
                               NodalGrad<dim> &, NodalGrad<dimLS> &,
                               SVec<double,3> &, Vec<double> &,
                               SVec<double,dim> &, GenMat<Scalar,neq> &,
                               FluidSelector &, Vec<int> &, int * );
  
  template<class Scalar,int dim,int neq>
  void computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>&,
                              FluxFcn**, GeoState&, SVec<double,3>&,
                              SVec<double,dim>&, Vec<double>&, LevelSetStructure &,
                              Vec<int>&, int, SVec<double,3>*,
                              GenMat<Scalar,neq>&,Vec<double>& irey);

  template<class Scalar,int dim, int dimLS,int neq>
  void computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,int*,
                                     FluxFcn** fluxFcn, 
                                     GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V, 
                                     LevelSetStructure& LSS, Vec<int> &fluidId,
                                     int Nriemann, SVec<double,3>* Nsbar, FluidSelector &fluidSelector,
                                     NodalGrad<dimLS>& ngradLS,Vec<double>&,GenMat<Scalar,neq>& A);

  template<class Scalar, int dim, int dimLS>
    void computeJacobianFiniteVolumeTermLS(RecFcn* recFcn, RecFcn* recFcnLS,
					   GeoState& geoState, SVec<double,3>& X,
					   SVec<double,dim>& V, NodalGrad<dim>& ngrad,
					   NodalGrad<dimLS> &ngradLS,
					   EdgeGrad<dim>* egrad,
					   Vec<double> &ctrlVol, SVec<double,dimLS>& Phi,
					   GenMat<Scalar,dimLS> &A,LevelSetStructure* LSS);
/*  template<int dim>
  void RiemannJacobianGasTait(int i, int j,
                              SVec<double,dim> &V, double Phii, double Phij,
		              double *nphi,
                              double *normal, double normalVel, VarFcn *varFcn,
                              FluxFcn** fluxFcn, double *dfdUi, double *dfdUj);
*/

  template<int dimLS>
  void TagInterfaceNodes(int lsdim, Vec<int> &Tag, SVec<double,dimLS> &Phi);
  template<int dimLS>
  void TagInterfaceNodes(int lsdim, Vec<int> &Tag1, Vec<int> &Tag2, SVec<double,dimLS> &Phi, LevelSetStructure *LSS);

  void setMasterFlag(bool *flag) { masterFlag = flag; }
  bool *getMasterFlag() const { return masterFlag; }
  int (*getPtr() const)[2] { return ptr; }
  int size() const { return numEdges; }

  void updateLength(SVec<double,3>& X);
  double length(int iedge) {
    if(edgeLength && iedge<numEdges) return(edgeLength[iedge]);
    else return(0.0);
  }
  double* viewEdgeLength() { return(edgeLength); }

// Included (MB)
  template<int dim>
  void computeDerivativeOfFiniteVolumeTerm(Vec<double> &, Vec<double> &, FluxFcn**, RecFcn*, ElemSet&, GeoState&, SVec<double,3>&, SVec<double,3>&,
			       SVec<double,dim>&, SVec<double,dim>&, NodalGrad<dim>&, EdgeGrad<dim>*, double,
			       SVec<double,dim>&);
  template<int dim>
  void computeDerivativeOfTimeStep(FemEquationTerm *, VarFcn *, GeoState &,
                              SVec<double,3> &, SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &,
			      Vec<double> &, Vec<double> &, double, TimeLowMachPrec &);
  int checkReconstructedValues(int i, int j, double *Vi, double *Vj, VarFcn *vf,
			       int *locToGlobNodeMap, int failsafe, SVec<int,2> &tag,
                               double *originalVi = 0, double *originalVj = 0,
                               int IDi = 0, int IDj = 0);
  void computeCharacteristicEdgeLength(SVec<double,3> &, double&, double&, double&, int&,
                                       const double, const double, const double,
                                       const double, const double, const double);

  void computeCellAveragedStructNormal(SVec<double,3> &, Vec<double> &, LevelSetStructure &);

  void attachProgrammedBurn(ProgrammedBurn*);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <Edge.C>
#endif

#endif
