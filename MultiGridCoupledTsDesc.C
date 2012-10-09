#include <VectorSet.h>
#include <LevelSet/LevelSetStructure.h>
#include <MultiGridCoupledTsDesc.h>

template <int dim>
MultiGridCoupledTsDesc<dim>::
MultiGridCoupledTsDesc(IoData & iod, GeoSource & gs,  Domain * dom) :
  ImplicitCoupledTsDesc<dim>(iod,gs,dom) {

}

template <int dim>
MultiGridCoupledTsDesc<dim>::
~MultiGridCoupledTsDesc() {

}

template class MultiGridCoupledTsDesc<5>;
template class MultiGridCoupledTsDesc<6>;
template class MultiGridCoupledTsDesc<7>;

