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

// Included
  SVec<Scalar,dim> *dddx;
  SVec<Scalar,dim> *dddy;
  SVec<Scalar,dim> *dddz;

public:

// Included
  NodalGrad(SVec<Scalar,dim> &dx, SVec<Scalar,dim> &dy, SVec<Scalar,dim> &dz,
            SVec<Scalar,dim> &dX, SVec<Scalar,dim> &dY, SVec<Scalar,dim> &dZ) :
            ddx(dx), ddy(dy), ddz(dz)
            { dddx = &dX; dddy = &dY; dddz = &dZ; }

  NodalGrad(SVec<Scalar,dim> &dx, SVec<Scalar,dim> &dy, SVec<Scalar,dim> &dz) : 
    ddx(dx), ddy(dy), ddz(dz) {}
  ~NodalGrad() {}

  SVec<Scalar,dim> &getX() const { return ddx; }
  SVec<Scalar,dim> &getY() const { return ddy; }
  SVec<Scalar,dim> &getZ() const { return ddz; }

// Included
  SVec<Scalar,dim> &getXderivative() const { return *dddx; }
  SVec<Scalar,dim> &getYderivative() const { return *dddy; }
  SVec<Scalar,dim> &getZderivative() const { return *dddz; }

};

//------------------------------------------------------------------------------

#endif
