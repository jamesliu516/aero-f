#include <Edge.h>

#include <Vector.h>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::swap;
#endif

#ifdef EDGE_LENGTH
#include <Vector3D.h>
#endif

//------------------------------------------------------------------------------

EdgeSet::EdgeSet()
{
  numEdges  = 0;
  ptr       = 0;
  masterFlag= 0;
  mp        = new MapPair;

#ifdef EDGE_LENGTH  //HB
  edgeLength= 0;
#endif
}

//------------------------------------------------------------------------------

EdgeSet::~EdgeSet()
{
  if (ptr)        delete [] ptr;
  if (masterFlag) delete [] masterFlag;
  if (mp)         delete mp;
#ifdef EDGE_LENGTH
  if (edgeLength) delete [] edgeLength;
#endif
}

//------------------------------------------------------------------------------

int EdgeSet::find(int first, int second)
{

  if (first > second) swap(first, second);

  MapPair::iterator v = mp->find(Pair(first, second));

  if (v == mp->end()) { 
    (*mp)[ Pair(first, second) ] = numEdges;
    return numEdges++;
  }

  return (*v).second;

}

//------------------------------------------------------------------------------

void EdgeSet::createPointers(Vec<int> &newNum)
{

  ptr = new int[numEdges][2];

  MapPair::iterator it = mp->begin();
  MapPair::iterator last = mp->end();

  int l = 0;

  while (it != last) {

    // XML I am not sure why one would need to renumber...
    ptr[l][0] = (*it).first.first;
    ptr[l][1] = (*it).first.second;

    newNum[(*it).second] = l;

    (*it).second = l;

    ++l;
    ++it;

  }

}

//------------------------------------------------------------------------------

#ifdef EDGE_LENGTH //HB: compute the edges'length for given nodes'coordinates
void EdgeSet::updateLength(SVec<double,3>& X)
{
  if(!edgeLength) edgeLength = new double[numEdges];
 
  for(int iedge=0; iedge<numEdges; iedge++) {
    Vec3D d(X[ptr[iedge][0]]);
    Vec3D e(X[ptr[iedge][1]]);
    Vec3D ed = d-e;
    edgeLength[iedge] = ed.norm();
  }
}
#endif

//------------------------------------------------------------------------------
int EdgeSet::checkReconstructedValues(int i, int j, double *Vi, double *Vj, VarFcn *vf, 
				int *locToGlobNodeMap, int failsafe, SVec<int,2>& tag,
                                double *originalVi, double *originalVj,
                                double phii, double phij)
{

  double rho = 0.0;
  double p   = 0.0;

// at interface of two-phase flow simulations, reverts to original value of pressure/density if
// the reconstructed ones are negative. 
// checkPressure does not check that pressure > 0, it checks that c^2>0, ie 
//     for stiffened gas, check that P+P_\infty>0
//     for Tait         , check not really needed but checks that P>0
/*  if(phii*phij<0){
    if (vf->getDensity(Vi)          <=0.0) vf->setDensity(Vi,originalVi);
    if (vf->checkPressure(Vi,phii) <= 0.0) vf->setPressure(Vi,originalVi,phii);

    if (vf->getDensity(Vj)          <=0.0) vf->setDensity(Vj,originalVj);
    if (vf->checkPressure(Vj,phij) <= 0.0) vf->setPressure(Vj,originalVj,phij);
  }
*/
//proceed to checking positivity of pressure and density for both nodes of an edge.
  int ierr = 0;

  rho = vf->getDensity(Vi);
  p   = vf->checkPressure(Vi,phii);

  if (rho <= 0.0) {
    if(!failsafe){
      fprintf(stderr, "*** Error: negative density (%e) for node %d after reconstruction on edge %d(%e) -> %d(%e)\n",
          rho, locToGlobNodeMap[i]+1, locToGlobNodeMap[i]+1, phii, locToGlobNodeMap[j]+1, phij);
      ++ierr;
    }
    else {
      fprintf(stderr, "*** Warning: negative density (%e) for node %d after reconstruction on edge %d(%e) -> %d(%e)\n",
           rho, locToGlobNodeMap[i]+1, locToGlobNodeMap[i]+1, phii, locToGlobNodeMap[j]+1, phij);
      tag[i][0] = 1;
      ++ierr;
    }
  }
  if (p <= 0.0) {
    if(!failsafe) {
      fprintf(stderr, "*** Error: negative pressure (%e) for node %d (rho = %e) after reconstruction on edge %d(%e) -> %d(%e)\n",
             p, locToGlobNodeMap[i]+1 , rho, locToGlobNodeMap[i]+1, phii, locToGlobNodeMap[j]+1, phij);
      ++ierr;
    }
    else {
     fprintf(stderr, "*** Warning: negative pressure (%e) for node %d (rho = %e) after reconstruction on edge %d(%e) -> %d(%e)\n",
             p, locToGlobNodeMap[i]+1 , rho, locToGlobNodeMap[i]+1, phii, locToGlobNodeMap[j]+1, phij);
     tag[i][0] = 1;
      ++ierr;
    }
  }



  rho = vf->getDensity(Vj);
  p   = vf->checkPressure(Vj,phij);

  if (rho <= 0.0) {
    if(!failsafe){
      fprintf(stderr, "*** Error: negative density (%e) for node %d after reconstruction on edge %d(%e) -> %d(%e)\n",
          rho, locToGlobNodeMap[j]+1, locToGlobNodeMap[j]+1, phij, locToGlobNodeMap[i]+1, phii);
      ++ierr;
    }
    else {
      fprintf(stderr, "*** Warning: negative density (%e) for node %d after reconstruction on edge %d(%e) -> %d(%e)\n",
          rho, locToGlobNodeMap[j]+1, locToGlobNodeMap[j]+1, phij, locToGlobNodeMap[i]+1, phii);
      tag[j][0] = 1;
      ++ierr;
    }
  }
  if (p <= 0.0) {
    if(!failsafe) {
      fprintf(stderr, "*** Error: negative pressure (%e) for node %d (rho = %e) after reconstruction on edge %d(%e) -> %d(%e)\n",
             p, locToGlobNodeMap[j]+1 , rho, locToGlobNodeMap[j]+1, phij, locToGlobNodeMap[i]+1, phii);
      ++ierr;
    }
    else {
     fprintf(stderr, "*** Warning: negative pressure (%e) for node %d (rho = %e) after reconstruction on edge %d(%e) -> %d(%e)\n",
             p, locToGlobNodeMap[j]+1 , rho, locToGlobNodeMap[j]+1, phij, locToGlobNodeMap[i]+1, phii);
     tag[j][0] = 1;
      ++ierr;
    }
  }



  return ierr;

}
//------------------------------------------------------------------------------
void EdgeSet::TagInterfaceNodes(Vec<int> &Tag, Vec<double> &Phi)
{

  int tag = 1;
  for(int l=0; l<numEdges; l++){
    int i = ptr[l][0];
    int j = ptr[l][1];

    if(Phi[i]*Phi[j]<=1.0e-12){
      Tag[i] = tag;
      Tag[j] = tag;
    }
  }

}
//------------------------------------------------------------------------------
