/*
 * DynamicNodalTransfer.h
 *
 *  Created on: May 12, 2009
 *      Author: Michel Lesoinne
 */

#ifndef DYNAMICNODALTRANSFER_H_
#define DYNAMICNODALTRANSFER_H_

class Communicator;
template<typename Scalar, int dim> class SVec;

/** Class to handle communication of nodal forces and displacement with the structure
 *
 */
class DynamicNodalTransfer {
	Communicator &com;
public:
	DynamicNodalTransfer(Communicator &);
	~DynamicNodalTransfer();

	void sendForce(SVec<double, 3> &);
};

#endif /* DYNAMICNODALTRANSFER_H_ */
