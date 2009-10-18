/*
 * DynamicNodalTransfer.h
 *
 *  Created on: May 12, 2009
 *      Author: Michel Lesoinne
 */

#ifndef DYNAMICNODALTRANSFER_H_
#define DYNAMICNODALTRANSFER_H_
#include<Vector.h>
#include<Communicator.h>
#include <DistInfo.h>
#include <map>

using std::pair;

class IoData;
class StructExc;

/** Class to temporarily play the role of a structure codes.
 *
 */
class EmbeddedStructure {
  Communicator &com;
  StructExc* structExc;

  const char *meshFile;
  const char *matcherFile;
 
  int mode;
  int coupled;
  double tScale;
  double dt, tMax;
  double omega;
  double dx, dy, dz;
  double t0; // starting time.
  int it;

  int nNodes;
  double (*X)[3]; //original node coordinates
  double (*U)[3]; //displacement
  double (*F)[3]; //force (received from fluid).
  std::map<int,int> pairing;

  DistInfo *di;
  
public:
  EmbeddedStructure(IoData& iod, Communicator &fc, Communicator &sc);
  ~EmbeddedStructure();

  pair<double*, int> getTargetData();
  void sendTimeStep(Communication::Window<double> *window);
  void sendDisplacement(Communication::Window<double> *window);
  void processReceivedForce();
};


/** Class to handle communication of nodal forces and displacement with the structure
 *
 */
class DynamicNodalTransfer {
        const double fScale; //reference force
        const double XScale; //reference length
        const double tScale; //reference time

	Communicator &com;
        EmbeddedStructure structure;

        SVec<double,3> F; //TODO: need to be resit by resetOutputToStructure
        double dts;
public:
	DynamicNodalTransfer(IoData& iod, Communicator &, Communicator &);
	~DynamicNodalTransfer();

	/** routine to send the force to the structure. On output, f has been dimensionalized. */
	void sendForce();
        /** routine to receive the displacement of the structure.*/
	void getDisplacement(SVec<double,3>& structU);

        double getStructureTimeStep() {return dts;}

        void updateOutputToStructure(double  dt, double dtLeft, SVec<double,3> &Fs);
};

#endif /* DYNAMICNODALTRANSFER_H_ */
