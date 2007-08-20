#ifndef _BC_DATA_H_
#define _BC_DATA_H_

#include <Vector.h>

//------------------------------------------------------------------------------

template<int dim>
class BcData {

  SVec<double,dim> &Uface;
  SVec<double,dim> &Unode;
  SVec<double,dim> &Uinletnode;
  SVec<double,dim> &Ufarin;
  SVec<double,dim> &Ufarout;

// Included (MB)
  SVec<double,dim> *dUface;
  SVec<double,dim> *dUnode;
  SVec<double,dim> *dUinletnode;
  SVec<double,dim> *dUfarin;
  SVec<double,dim> *dUfarout;
  SVec<double,dim> *dUfaceSA;
  SVec<double,dim> *dUnodeSA;

public:

// Included (MB)
  BcData(SVec<double,dim> &uf, SVec<double,dim> &un, SVec<double,dim> &uin, SVec<double,dim> &ufarin, SVec<double,dim> &ufarout, SVec<double,dim> &duf, SVec<double,dim> &dun, SVec<double,dim> &duin, SVec<double,dim> &dufarin, SVec<double,dim> &dufarout, SVec<double,dim> &dufsa, SVec<double,dim> &dunsa)
         : Uface(uf), Unode(un), Uinletnode(uin), Ufarin(ufarin), Ufarout(ufarout) { dUface = &duf; dUnode = &dun; dUinletnode = &duin; dUfarin = &dufarin; dUfarout = &dufarout; dUfaceSA = &dufsa; dUnodeSA = &dunsa;}
  BcData(SVec<double,dim> &uf, SVec<double,dim> &un, SVec<double,dim> &uin, SVec<double,dim> &ufarin, SVec<double,dim> &ufarout, SVec<double,dim> &duf, SVec<double,dim> &dun, SVec<double,dim> &duin, SVec<double,dim> &dufarin, SVec<double,dim> &dufarout)
         : Uface(uf), Unode(un), Uinletnode(uin), Ufarin(ufarin), Ufarout(ufarout) { dUface = &duf; dUnode = &dun; dUinletnode = &duin; dUfarin = &dufarin; dUfarout = &dufarout; }
  BcData(SVec<double,dim> &uf, SVec<double,dim> &un, SVec<double,dim> &uin, SVec<double,dim> &ufarin, SVec<double,dim> &ufarout, SVec<double,dim> &dufsa, SVec<double,dim> &dunsa)
         : Uface(uf), Unode(un), Uinletnode(uin), Ufarin(ufarin), Ufarout(ufarout) { dUfaceSA = &dufsa; dUnodeSA = &dunsa; }

  BcData(SVec<double,dim> &uf, SVec<double,dim> &un, SVec<double,dim> &uin, SVec<double,dim> &ufarin, SVec<double,dim> &ufarout) : Uface(uf), Unode(un), Uinletnode(uin), Ufarin(ufarin), Ufarout(ufarout) {}
  ~BcData() {}

  SVec<double,dim> &getFaceStateVector() const { return Uface; }
  SVec<double,dim> &getNodeStateVector() const { return Unode; }
  SVec<double,dim> &getInletNodeStateVector() const { return Uinletnode; }
  SVec<double,dim> &getInletBoundaryVector() const { return Ufarin; }
  SVec<double,dim> &getInletOutletVector() const { return Ufarout; }

// Included (MB)
  SVec<double,dim> &getdFaceStateVector() const { return *dUface; }
  SVec<double,dim> &getdNodeStateVector() const { return *dUnode; }
  SVec<double,dim> &getdFaceStateVectorSA() const { return *dUfaceSA; }
  SVec<double,dim> &getdNodeStateVectorSA() const { return *dUnodeSA; }
  SVec<double,dim> &getdInletNodeStateVector() const { return *dUinletnode; }
  SVec<double,dim> &getdInletBoundaryVector() const { return *dUfarin; }
  SVec<double,dim> &getdInletOutletVector() const { return *dUfarout; }

};

//------------------------------------------------------------------------------

#endif
