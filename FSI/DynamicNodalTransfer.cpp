/*
 * DynamicNodalTransfer.cpp
 *
 *  Created on: May 12, 2009
 *      Author: michel
 */
#include <iostream>
#include <FSI/DynamicNodalTransfer.h>
#include <Vector.h>
#include <Communicator.h>
#include <IoData.h>

DynamicNodalTransfer::DynamicNodalTransfer(IoData& iod, Communicator &c): com(c) {
  fScale = iod.ref.rv.tforce;
  std::cout << "fScale is " << fScale << std::endl;
}

DynamicNodalTransfer::~DynamicNodalTransfer() {

}

void
DynamicNodalTransfer::sendForce(SVec<double, 3> &f) {
  f *= fScale;
  std::cout << "Sending a force with a norm " << f.norm() << std::endl;
  Communication::Window<double> window(com, 3*f.size()*sizeof(double), (double *)f.data());
  window.fence(true);
  window.accumulate((double *)f.data(), 0, 3*f.size(), 0, 0, Communication::Window<double>::Add);
  window.fence(false);
}
