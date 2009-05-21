/*
 * DynamicNodalTransfer.h
 *
 *  Created on: May 12, 2009
 *      Author: Michel Lesoinne
 */

#ifndef DYNAMICNODALTRANSFER_H_
#define DYNAMICNODALTRANSFER_H_

class Communicator;
class IoData;
template<typename Scalar, int dim> class SVec;

/** Class to handle communication of nodal forces and displacement with the structure
 *
 */
class DynamicNodalTransfer {
  double fScale;
	Communicator &com;
public:
	DynamicNodalTransfer(IoData& iod, Communicator &);
	~DynamicNodalTransfer();

	/** routine to send the force to the structure. On output, f has been dimensionalized. */
	void sendForce(SVec<double, 3> &f);
};

#endif /* DYNAMICNODALTRANSFER_H_ */
