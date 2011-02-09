#ifndef _STRUCT_EXC_H_
#define _STRUCT_EXC_H_

class IoData;
class GeoSource;
class MatchNodeSet;
class Domain;
class Communicator;

template<class Scalar> class DistVec;
template<class Scalar, int dim> class DistSVec;

//------------------------------------------------------------------------------

class StructExc {

  int algNum;
  int rstrt;
  int smode;
  double dt;
  double tmax;

  double oolscale;
  double ootscale;
  double oovscale;
  double ootempscale;
  double fscale;
  double pscale;

  int numLocSub;
  int numStrCPU;
  int (*numStrNodes)[2];

  int recParity;
  int sndParity;
  int bufsize;
  double *buffer;

  MatchNodeSet **matchNodes;

  Communicator *com;
  Communicator *strCom;

public:

  StructExc(IoData&, MatchNodeSet**, int, Communicator*strCom, Communicator *flCom, int nSub);
  ~StructExc();
  void updateMNS(MatchNodeSet **mns) {matchNodes = mns;}
  void updateNumStrNodes(int nn);

  void negotiate();
  void negotiateStopping(bool*);
  double getInfo();
  void getEmbeddedWetSurfaceInfo(int&, bool&, int&, int&);
  void getEmbeddedWetSurface(int, double*, int, int*, int=3);  //3 for triangle

  void getDisplacement(DistSVec<double,3> &, DistSVec<double,3> &, 
		       DistSVec<double,3> &, DistSVec<double,3> &);
  void getTemperature(DistVec<double>&);
  void sendForce(DistSVec<double,3> &);
  void getMdFreq(int &, double *&);
  void getMdStrDisp(int, DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,3> &);
  void sendHeatPower(DistVec<double>&);

  int getAlgorithmNumber() const { return algNum; }
  int getRestartFrequency() const { return rstrt; }
  double getTimeStep() const { return dt; }
  double getMaxTime() const { return tmax; }

  //for cracking
  void getInitialCrackingSetup(int&, int&);
  int getNumberOfNewCrackedElems();
  void getInitialPhantomNodes(int nCracked, double(*xyz)[3], int nNodes);
  void getNewCracking(int nCracked, int* phantoms, double* phi);
};

//------------------------------------------------------------------------------

#endif
