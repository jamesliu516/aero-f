#pragma once

#include <ImplicitCoupledTsDesc.h>

template <int dim>
class MultiGridCoupledTsDesc : public ImplicitCoupledTsDesc<dim> {

 public:

  MultiGridCoupledTsDesc(IoData &, GeoSource &, Domain *); 

  ~MultiGridCoupledTsDesc();
};
