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

class VarFcn;
class RecFcn;
class FluxFcn;
class TetSet;
class GeoState;
class FemEquationTerm;

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

#ifdef EDGE_LENGTH
#include <Vector.h>
#endif

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

#ifdef EDGE_LENGTH  //HB: stores the (current) length of the edges
  double* edgeLength;
#endif

public:

  EdgeSet();
  ~EdgeSet();

  int find(int, int);

  void createPointers(Vec<int> &);

  template<int dim>
  void computeTimeStep(FemEquationTerm *, VarFcn *, GeoState &, 
			SVec<double,3> &, SVec<double,dim> &, Vec<double> &,
                       Vec<double> &, double, double, double);
  template<int dim>
  void computeTimeStep(VarFcn *, GeoState &, SVec<double,dim> &, Vec<double> &,
                       double, double, double, Vec<double> &, int);

  template<int dim>
  int computeFiniteVolumeTerm(int*, Vec<double> &, FluxFcn**, RecFcn*, TetSet&, GeoState&, 
                              SVec<double,3>&, SVec<double,dim>&, NodalGrad<dim>&, EdgeGrad<dim>*,
			      SVec<double,dim>&, SVec<int,2>&, int, int);

  template<int dim>
  int computeFiniteVolumeTerm(ExactRiemannSolver<dim>&, int*,
                              FluxFcn**, RecFcn*, TetSet&, GeoState&, SVec<double,3>&,
                              SVec<double,dim>&, Vec<double> &, 
                              NodalGrad<dim>&, EdgeGrad<dim>*,
                              NodalGrad<1>&,
                              SVec<double,dim>&, int,
                              SVec<int,2>&, int, int);

  template<int dim>
  void computeFiniteVolumeTermLS(FluxFcn**, RecFcn*, RecFcn*, TetSet&, GeoState&, SVec<double,3>&,
                               SVec<double,dim>&, NodalGrad<dim>&, NodalGrad<1>&, EdgeGrad<dim>*,
                               SVec<double,1>&, Vec<double>&);

  template<int dim>
  void computeFiniteVolumeBarTerm(FluxFcn**, RecFcn*, TetSet&, GeoState&, SVec<double,3>&,
			       SVec<double,dim>&, NodalGrad<dim>&, EdgeGrad<dim>*,
			       SVec<double,dim>&, SVec<double,1>&);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, GeoState &, 
				Vec<double> &, Vec<double> &,
			       SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, GeoState &,
				Vec<double> &, Vec<double> &,
                               SVec<double,dim> &, GenMat<Scalar,neq> &, 
                               int * );

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, GeoState &, 
 			       NodalGrad<dim> &, NodalGrad<1> &, Vec<double> &,
			       SVec<double,dim> &, GenMat<Scalar,neq> &,
                               Vec<double> &);
  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, GeoState &, 
 			       NodalGrad<dim> &, NodalGrad<1> &, Vec<double> &,
                               SVec<double,dim> &, GenMat<Scalar,neq> &,
                               Vec<double> &, int * );

/*  template<int dim>
  void RiemannJacobianGasTait(int i, int j,
                              SVec<double,dim> &V, double Phii, double Phij,
		              double *nphi, 
                              double *normal, double normalVel, VarFcn *varFcn,
                              FluxFcn** fluxFcn, double *dfdUi, double *dfdUj);
*/

  void TagInterfaceNodes(Vec<int> &Tag, Vec<double> &Phi);

  void setMasterFlag(bool *flag) { masterFlag = flag; }
  bool *getMasterFlag() const { return masterFlag; }
  int (*getPtr() const)[2] { return ptr; }
  int size() const { return numEdges; }

#ifdef EDGE_LENGTH  
  void updateLength(SVec<double,3>& X);
  double length(int iedge) { if(edgeLength && iedge<numEdges) return(edgeLength[iedge]);
                             else return(0.0); }
  double* viewEdgeLength() { return(edgeLength); }
#endif

  int checkReconstructedValues(int i, int j, double *Vi, double *Vj, VarFcn *vf,
			       int *locToGlobNodeMap, int failsafe, SVec<int,2> &tag,
                               double phii = 1.0, double phij = 1.0);

  int checkReconstruction(double [2], double [2], int, int, int*, int, SVec<int,2>&);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <Edge.C>
#endif

#endif
