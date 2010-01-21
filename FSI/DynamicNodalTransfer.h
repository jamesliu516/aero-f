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
#include <Timer.h>
#include <map>

using std::pair;

class IoData;
class StructExc;

/** Class to temporarily play the role of a structure codes.
 *
 */
class EmbeddedStructure {
  Communicator &com;
  Timer *timer;
  StructExc* structExc;
  
  const char *meshFile;
  const char *matcherFile;
 
  bool coupled;
  int algNum;
  bool dim2Treatment;
  bool oneWayCoupling;
  int mode;
  
  double tScale;
  double XScale;
  double UScale;

  double dt, tMax;
  double omega;
  double dx, dy, dz;
  double t0; // starting time.
  int it;

  int nNodes;
  double (*X)[3]; //original node coordinates
  double (*U)[3]; //displacement
  double (*Udot)[3]; //velocity
  double (*UandUdot)[3]; //displacement and velocity (TODO: this is redundant)
  double (*F)[3]; //force (received from fluid).
  std::map<int,int> pairing;

  DistInfo *di;
  
public:
  EmbeddedStructure(IoData& iod, Communicator &fc, Communicator &sc, Timer *tim);
  ~EmbeddedStructure();

  int getAlgorithmNumber() {return algNum;}
  int numStructNodes() {return nNodes;}
  pair<double*, int> getTargetData();
  void sendTimeStep(Communication::Window<double> *window);
  void sendMaxTime(Communication::Window<double> *window);
  void sendDisplacement(Communication::Window<double> *window);
  void processReceivedForce();
};


/** Class to handle communication of nodal forces and displacement with the structure
 *
 */
class DynamicNodalTransfer {
        const double fScale; //scaling factor for force
        const double XScale; //scaling factor for length
        const double tScale; //scaling factor for time
        const double UScale; //scaling factor for velocity
        int algNum;

	Communicator &com;
        Timer *timer;
        EmbeddedStructure structure;

        SVec<double,3> F; //TODO: need to be resit by resetOutputToStructure
        double dts;
        double tMax;
public:
	DynamicNodalTransfer(IoData& iod, Communicator &, Communicator &, Timer *);
	~DynamicNodalTransfer();

        int getAlgorithmNumber() {return algNum;}

	/** routine to send the force to the structure. On output, f has been dimensionalized. */
	void sendForce();
        /** routine to receive the displacement of the structure.*/
	void getDisplacement(SVec<double,3>& structU, SVec<double,3>& structUdot);

        double getStructureTimeStep() {return dts;}
        double getStructureMaxTime() {return tMax;}

        void updateOutputToStructure(double  dt, double dtLeft, SVec<double,3> &Fs);
};

#endif /* DYNAMICNODALTRANSFER_H_ */
