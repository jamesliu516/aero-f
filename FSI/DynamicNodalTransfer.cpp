/*
 * DynamicNodalTransfer.cpp
 *
 *  Created on: May 12, 2009
 *      Author: michel
 */
#include <iostream>
#include <FSI/DynamicNodalTransfer.h>
#include <Communicator.h>
#include <IoData.h>

DynamicNodalTransfer::DynamicNodalTransfer(IoData& iod, Communicator &c): com(c) , F(1), 
                           fScale(iod.ref.rv.tforce), XScale(iod.ref.rv.tlength),
                           tScale(iod.ref.rv.time)

{
  fprintf(stderr,"fscale = %e, XScale = %e, tScale = %e.\n", fScale, XScale, tScale);
  Communication::Window<double> window(com, 1, &dts);
  window.fence(true);
  window.fence(false);
  dts /= tScale;
  fprintf(stderr,"dt = %e\n", dts);

}

DynamicNodalTransfer::~DynamicNodalTransfer() {

}

void
DynamicNodalTransfer::sendForce() {
  Communication::Window<double> window(com, 3*F.size()*sizeof(double), (double *)F.data());
  window.fence(true);
  window.accumulate((double *)F.data(), 0, 3*F.size(), 0, 0, Communication::Window<double>::Add);
  window.fence(false);
}

void
DynamicNodalTransfer::getDisplacement(SVec<double,3>& structU) {
  Communication::Window<double> window(com, 3*structU.size()*sizeof(double), (double *)structU.data());
  window.fence(true);
  window.fence(false);

  structU = 1.0/XScale*structU;
}

double
DynamicNodalTransfer::getStructureTimeStep() {

 return dts;
} 

void
DynamicNodalTransfer::updateOutputToStructure(double dt, double dtLeft, SVec<double,3> &fs)
{
  if(F.size() != fs.size())
    F.resize(fs.size());
  F = fScale * fs;
}
