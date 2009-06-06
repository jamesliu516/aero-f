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

  void negotiate();
  void negotiateStopping(bool*);
  double getInfo();
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

};

//------------------------------------------------------------------------------

#endif
