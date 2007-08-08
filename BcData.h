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

public:

  BcData(SVec<double,dim> &uf, SVec<double,dim> &un, SVec<double,dim> &uin, SVec<double,dim> &ufarin, SVec<double,dim> &ufarout) : Uface(uf), Unode(un), Uinletnode(uin), Ufarin(ufarin), Ufarout(ufarout) {}
  ~BcData() {}

  SVec<double,dim> &getFaceStateVector() const { return Uface; }
  SVec<double,dim> &getNodeStateVector() const { return Unode; }
  SVec<double,dim> &getInletNodeStateVector() const { return Uinletnode; }
  SVec<double,dim> &getInletBoundaryVector() const { return Ufarin; }
  SVec<double,dim> &getInletOutletVector() const { return Ufarout; }

};

//------------------------------------------------------------------------------

#endif
