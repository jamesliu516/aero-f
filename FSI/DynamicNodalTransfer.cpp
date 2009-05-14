/*
 * DynamicNodalTransfer.cpp
 *
 *  Created on: May 12, 2009
 *      Author: michel
 */

#include <FSI/DynamicNodalTransfer.h>
#include <Vector.h>
#include <Communicator.h>

DynamicNodalTransfer::DynamicNodalTransfer(Communicator &c): com(c) {

}

DynamicNodalTransfer::~DynamicNodalTransfer() {

}

void
DynamicNodalTransfer::sendForce(SVec<double, 3> &f) {
  Communication::Window<double> window(com, 3*f.size()*sizeof(double), (double *)f.data());
  window.fence(true);
  window.accumulate((double *)f.data(), 0, 3*f.size(), 0, 0, Communication::Window<double>::Add);
  window.fence(false);
}
