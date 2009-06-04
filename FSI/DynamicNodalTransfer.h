/*
 * DynamicNodalTransfer.h
 *
 *  Created on: May 12, 2009
 *      Author: Michel Lesoinne
 */

#ifndef DYNAMICNODALTRANSFER_H_
#define DYNAMICNODALTRANSFER_H_
#include<Vector.h>

class Communicator;
class IoData;

/** Class to handle communication of nodal forces and displacement with the structure
 *
 */
class DynamicNodalTransfer {
        const double fScale; //reference force
        const double XScale; //reference length
        const double tScale; //reference time

	Communicator &com;
        SVec<double,3> F; //TODO: need to be resit by resetOutputToStructure
        double dts;
public:
	DynamicNodalTransfer(IoData& iod, Communicator &);
	~DynamicNodalTransfer();

	/** routine to send the force to the structure. On output, f has been dimensionalized. */
	void sendForce();
        /** routine to receive the displacement of the structure.*/
	void getDisplacement(SVec<double,3>& structU);

        double getStructureTimeStep();

        void updateOutputToStructure(double  dt, double dtLeft, SVec<double,3> &Fs);
};

#endif /* DYNAMICNODALTRANSFER_H_ */
