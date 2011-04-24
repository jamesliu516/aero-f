/*
 * DynamicNodalTransfer.h
 *
 *  Created on: May 12, 2009
 *      Author: Michel Lesoinne, Kevin Wang
 */

#ifndef DYNAMICNODALTRANSFER_H_
#define DYNAMICNODALTRANSFER_H_
#include <Vector.h>
#include <Communicator.h>
#include <DistInfo.h>
#include <Timer.h>
#include <map>
#include <string.h>

using std::pair;

class IoData;
class StructExc;
class CrackingSurface;
class MatchNodeSet;

class EmbeddedStructure {
  friend class DynamicNodalTransfer;

  Communicator &com;
  Timer *timer;
  StructExc* structExc;
  
  bool getSurfFromFEM;

  char *meshFile;
  char *matcherFile;
 
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

  int nNodes, totalNodes;
  int nElems, totalElems;
  int elemType;
  CrackingSurface *cracking; //activated only if cracking is considered in the structure code.

  //NOTE: the following variables should be dimensional!
  double (*X0)[3]; //original node coordinates
  double (*X)[3];  //updated node coordinates
  int    (*Tria)[3]; //mesh topology (activated only if the mesh is provided by FEM)
  double (*U)[3]; //displacement
  double (*Udot)[3]; //velocity
  double (*XandUdot)[3]; //displacement and velocity (TODO: this is redundant)
  double dt_tmax[2];
  double (*F)[3]; //force (received from fluid).
  std::map<int,int> pairing;

  DistInfo *di;
  MatchNodeSet **mns;
  
  void splitQuads(int*,int); //utility function which split quads into triangles
  void getInitialCrack();
  int getNewCracking();
  void getNewCracking(int,int,int);
  void clearForceVector(); 
  pair<double*, int> getTargetData();

public:
  EmbeddedStructure(IoData& iod, Communicator &fc, Communicator &sc, Timer *tim);
  ~EmbeddedStructure();

  int getAlgorithmNumber() {return algNum;}
  int numStructNodes() {return nNodes;}
  void sendTimeStep(Communication::Window<double> *window);
  void sendMaxTime(Communication::Window<double> *window);
  void sendInfo(Communication::Window<double> *window);
  void sendInitialPosition(Communication::Window<double> *window);
  void sendDisplacement(Communication::Window<double> *window);
  void processReceivedForce();

  //if embedded mesh provided by FEM
  bool embeddedMeshByFEM() {return getSurfFromFEM;}
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

        Communication::Window<double> *wintime;
        Communication::Window<double> *winForce;
        Communication::Window<double> *winDisp;
        double *XandUdot; //non-dimensionalized
        double *dt_tmax;

        SVec<double,3> F; //TODO: need to be resit by resetOutputToStructure
        double dts;
        double tMax;
public:
	DynamicNodalTransfer(IoData& iod, Communicator &, Communicator &, Timer *);
	~DynamicNodalTransfer();

        int getAlgorithmNumber() {return algNum;}

	/** routine to send the force to the structure. On output, f has been dimensionalized. */
	void sendForce();
        /** routine to receive the new time step and final time from structure. */
        void updateInfo();
        /** routine to receive the new cracking of the structure.*/
	int getNewCracking();
        /** routine to receive the displacement of the structure.*/
	void getDisplacement();

        double getStructureTimeStep() {return dts;}
        double getStructureMaxTime() {return tMax;}

        void updateOutputToStructure(double  dt, double dtLeft, SVec<double,3> &Fs);
        bool embeddedMeshByFEM() {return structure.embeddedMeshByFEM();}
        int  numStNodes() {return structure.nNodes;}
        int  numStElems() {return structure.nElems;}
        int  totStNodes() {return structure.totalNodes;}
        int  totStElems() {return structure.totalElems;}
        bool cracking()   {return (structure.cracking) ? true : false;}
        double *getStNodes() {return XandUdot;}
        double *getStVelocity() {return XandUdot+structure.nNodes*3;}
        int    (*getStElems())[3] {return structure.Tria;}
        
        CrackingSurface* getCrackingSurface() {return structure.cracking;}

//        void getEmbeddedMesh(int &n1, double (**xyz)[3], int &n2, int (**abc)[3]) {
//          structure.getEmbeddedMesh(n1,xyz,n2,abc); 
//          fprintf(stderr,"DY %d %e %e %e\n", 2, xyz[1][0], xyz[1][1], xyz[1][2]);
//          fprintf(stderr,"DY %d %e %e %e\n", 2, structure.X[1][0], structure.X[1][1], structure.X[1][2]);}

};

#endif /* DYNAMICNODALTRANSFER_H_ */
