#ifndef _NODAL_GRAD_H_
#define _NODAL_GRAD_H_

#include <Vector.h>


#ifndef _NDGRAD_TMPL_
#define _NDGRAD_TMPL_
template<int dim, class Scalar = double> class NodalGrad;
#endif

//------------------------------------------------------------------------------

template<int dim, class Scalar>
class NodalGrad {

  SVec<Scalar,dim> &ddx;
  SVec<Scalar,dim> &ddy;
  SVec<Scalar,dim> &ddz;

public:

  NodalGrad(SVec<Scalar,dim> &dx, SVec<Scalar,dim> &dy, SVec<Scalar,dim> &dz) : 
    ddx(dx), ddy(dy), ddz(dz) {}
  ~NodalGrad() {}

  SVec<Scalar,dim> &getX() const { return ddx; }
  SVec<Scalar,dim> &getY() const { return ddy; }
  SVec<Scalar,dim> &getZ() const { return ddz; }

};

//------------------------------------------------------------------------------

#endif
