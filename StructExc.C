#include <StructExc.h>

#include <IoData.h>
#include <MatchNode.h>
#include <Domain.h>
#include <DistVector.h>

#include <stdlib.h>
#include <math.h>

#define FORCE_TAG 1000
#define DISP_TAG 2000
#define INFO_TAG 3000
#define HEATPOWER_TAG 5000
#define TEMP_TAG 6000
#define STRUC_CMD_TAG 8000
#define FLUID_CMD_TAG 9000
#define NEGO_NUM_TAG 10000
#define NEGO_BUF_TAG 10001

//------------------------------------------------------------------------------

StructExc::StructExc(IoData& iod, MatchNodeSet** mns, int bs, Communicator* sc, Domain* domain)
{

  bufsize = bs;
  com = domain->getCommunicator();
  strCom = sc;

  if (!strCom) {
    com->fprintf(stderr, "*** Error: structure communicator is null\n");
    exit(1);
  }

  oolscale = 1.0 / iod.ref.rv.tlength;
  ootscale = 1.0 / iod.ref.rv.time;
  oovscale = 1.0 / iod.ref.rv.tvelocity;
  ootempscale = 1.0 / iod.ref.rv.temperature;
  fscale = iod.ref.rv.tforce;
  pscale = iod.ref.rv.tpower;

  numLocSub = domain->getNumLocSub();
  numStrCPU = strCom->remoteSize();
  numStrNodes = 0;

  sndParity = 0;
  recParity = 0;
  buffer = 0;

  matchNodes = mns;

}

//------------------------------------------------------------------------------

StructExc::~StructExc()
{

  if (numStrNodes) delete [] numStrNodes;
  if (buffer) delete [] buffer;

}

//------------------------------------------------------------------------------

void StructExc::negotiate()
{

  // compute the total number of matched nodes for this fluid CPU

  int (*ptr)[2] = new int[numLocSub][2];

  int iSub;
#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub)
    ptr[iSub][0] = matchNodes[iSub]->size();

  ptr[0][1] = 0;
  for (iSub=1; iSub<numLocSub; ++iSub)
    ptr[iSub][1] = ptr[iSub - 1][1] + ptr[iSub - 1][0];

  int numCPUMatchedNodes = ptr[numLocSub - 1][1] + ptr[numLocSub - 1][0];

  // gather the global matched point nodes for this fluid CPU

  int *ibuffer = 0;
  int (*cpuMatchNodes)[3] = 0;

  if (numCPUMatchedNodes > 0) {
    ibuffer = new int[numCPUMatchedNodes];
    cpuMatchNodes = new int[numCPUMatchedNodes][3];
    buffer = new double[bufsize * numCPUMatchedNodes];

#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub)
      matchNodes[iSub]->exportInfo(iSub, cpuMatchNodes + ptr[iSub][1]);
  }

  // send the matched node numbers of this fluid CPU to all the structure CPUs

  int i;
  for (i=0; i<numCPUMatchedNodes; ++i)
    ibuffer[i] = cpuMatchNodes[i][1];

  int iCpu;
  for (iCpu=0; iCpu<numStrCPU; ++iCpu) {
    strCom->sendTo(iCpu, NEGO_NUM_TAG, &numCPUMatchedNodes, 1);
    if (numCPUMatchedNodes > 0)
      strCom->sendTo(iCpu, NEGO_BUF_TAG, ibuffer, numCPUMatchedNodes);
  }
  strCom->waitForAllReq();

  // receive the list of matched nodes that each structure CPU contains

  if (numCPUMatchedNodes > 0) {

    numStrNodes = new int[numStrCPU][2];

    int pos = 0;

    for (iCpu=0; iCpu<numStrCPU; ++iCpu) {
      strCom->recFrom(iCpu, NEGO_NUM_TAG, &numStrNodes[iCpu][0], 1);

      if (numStrNodes[iCpu][0] > 0) {
	strCom->recFrom(iCpu, NEGO_BUF_TAG, ibuffer, numStrNodes[iCpu][0]);
	for (i=0; i<numStrNodes[iCpu][0]; ++i) {
	  int idx = ibuffer[i];
	  matchNodes[cpuMatchNodes[idx][2]]->setBufferPosition(cpuMatchNodes[idx][0], pos);
	  pos++;
	}
      }

    }

    if (pos != numCPUMatchedNodes) {
      fprintf(stderr, "*** Error: wrong number of matched nodes (%d instead of %d)\n",
	      pos, numCPUMatchedNodes);
      exit(1);
    }

    numStrNodes[0][1] = 0;
    for (iCpu=1; iCpu<numStrCPU; ++iCpu)
      numStrNodes[iCpu][1] = numStrNodes[iCpu-1][1] + numStrNodes[iCpu-1][0];

  }

  if (com->getMaxVerbose() >= 8)
    fprintf(stdout, "CPU %d has %d matched node%s\n", 
	    com->cpuNum(), numCPUMatchedNodes, numCPUMatchedNodes>1? "s":"");

  if (ptr) delete [] ptr;
  if (ibuffer) delete [] ibuffer;
  if (cpuMatchNodes) delete [] cpuMatchNodes;

}

//------------------------------------------------------------------------------

double StructExc::getInfo() 
{

  double info[5];

  if (strCom->cpuNum() == 0)
    strCom->recFrom(INFO_TAG, info, 5);

  com->broadcast(5, info);

  algNum = int(info[0]);
  dt = info[1] * ootscale;
  tmax = info[2] * ootscale;
  rstrt = int(info[3]);
  smode = int(info[4]);

  if (algNum == 6) tmax -= 0.5 * dt;
  if (algNum == 20) tmax -= 0.5 * dt;
  if (algNum == 21) tmax += 1.5 * dt;

  double mppFactor = 1.0;
  if (algNum == 8)
    mppFactor = 1.0/info[2];

  return mppFactor;

}

//------------------------------------------------------------------------------
/*
     command communication between structure "0" and fluid "0"
     ! structure allways starts communication !

         structure sends            fluid has             fluid sends
     ------------------------------------------------------------------------
     small - talk

  1.     0 : continue               0 : continue          0: continue
  2.     0 : continue               1 : converged         0: continue
  3.     1 : converged              0 : continue          0: continue
  4.     1 : converged              1 : converged         1: stop analysis
                                                             next task
  5.     2 : converged              0/1                   1: stop analysis
                                                             next task
  5.    -1 : exit                   0/1                  -1: exit
  6.     0/1                       -1 : exit             -1: exit
     ------------------------------------------------------------------------
     command - talk

  7.   100 : analysis               0/1                 100: analysis
  8.   200 : senstivity             0/1                 200: senstivitiy
  9.   100/200                     -1 : exit            - 1: exit

*/

void StructExc::negotiateStopping(bool* lastIt)
{

  double xbuf;
  
  if (strCom->cpuNum() == 0) {
    strCom->recFrom(STRUC_CMD_TAG, &xbuf, 1);
    if (xbuf == 2.0)
      xbuf = 1.0;
    strCom->sendTo(0, FLUID_CMD_TAG, &xbuf, 1);
    strCom->waitForAllReq();
  }

  com->broadcast(1, &xbuf);  
  if (xbuf == 0.0)
    *lastIt = false;
  else
    *lastIt = true;

}

//------------------------------------------------------------------------------
/* 
   dX contains the displacement of the boundaries with respect to the CURRENT 
   configuration of the structure
*/

void StructExc::getDisplacement(DistSVec<double,3> &X0, DistSVec<double,3> &X, 
				DistSVec<double,3> &Xdot, DistSVec<double,3> &dX) 
{  

  double norms[2] = {0.0, 0.0};

  dX = 0.0;
  Xdot = 0.0;

  if (algNum == 4 || algNum == 5) recParity = 1 - recParity;

  if (numStrNodes) {

    for (int iCpu=0; iCpu<numStrCPU; ++iCpu) {
      if (numStrNodes[iCpu][0] > 0) {
	int size = bufsize * numStrNodes[iCpu][0];
	double *localBuffer = buffer + bufsize * numStrNodes[iCpu][1];
	strCom->recFrom(iCpu, DISP_TAG + recParity, localBuffer, size);
      }
    }

    double (*disp)[2][3] = reinterpret_cast<double (*)[2][3]>(buffer);

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      double (*x0)[3] = X0.subData(iSub);
      double (*x)[3] = X.subData(iSub);
      double (*xdot)[3] = Xdot.subData(iSub);
      double (*dx)[3] = dX.subData(iSub);

      double locNorms[2];
      matchNodes[iSub]->getDisplacement(algNum, dt, oolscale, oovscale, X.getMasterFlag(iSub), 
					disp, x0, x, xdot, dx, locNorms);
#pragma omp critical
      norms[0] += locNorms[0];
#pragma omp critical
      norms[1] += locNorms[1];
    }

  }

  com->globalSum(2, norms);

  com->printf(7, "Received total disp=%e and vel=%e from the structure\n", norms[0], norms[1]);
  com->printf(1, "Received total disp=%e and vel=%e from the structure\n", norms[0], norms[1]);
}

//------------------------------------------------------------------------------

void StructExc::getTemperature(DistVec<double>& Temp)
{  

  double norm = 0.0;

  if (numStrNodes) {
    for (int iCpu=0; iCpu<numStrCPU; ++iCpu) {
      if (numStrNodes[iCpu][0] > 0) {
	int size = bufsize * numStrNodes[iCpu][0];
	double* localBuffer = buffer + bufsize * numStrNodes[iCpu][1];
	strCom->recFrom(iCpu, TEMP_TAG, localBuffer, size);
      }
    }

#pragma omp parallel for reduction(+: norm)
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      double* temp = Temp.subData(iSub);
      bool* flag = Temp.getMasterFlag(iSub);
      norm += matchNodes[iSub]->getTemperature(algNum, dt, ootempscale, flag, buffer, temp);
    }
  }

  com->globalSum(1, &norm);

  com->printf(7, "Received temp=%e from the structure\n", norm);

}

//------------------------------------------------------------------------------
// note: the force vector is *** NOT *** assembled

void StructExc::sendForce(DistSVec<double,3> &F) 
{

  if (algNum == 4 || algNum == 5) sndParity = 1 - sndParity;

  double norm = 0.0;

  if (numStrNodes) {

    double (*forces)[3] = reinterpret_cast<double (*)[3]>(buffer);

#pragma omp parallel for reduction (+: norm)
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      SVec<double,3> &f = F(iSub);
      norm += f*f * fscale*fscale;
      matchNodes[iSub]->send(fscale, f.data(), forces);
    }

    for (int iCpu=0; iCpu<numStrCPU; ++iCpu) {
      if (numStrNodes[iCpu][0] > 0) {
	int size = 3 * numStrNodes[iCpu][0];
	double *localBuffer = buffer + 3 * numStrNodes[iCpu][1];
	strCom->sendTo(iCpu, FORCE_TAG + sndParity, localBuffer, size);
      }
    }
    strCom->waitForAllReq();

  }

  com->globalSum(1, &norm);

  com->printf(7, "Sent fluid force=%e to the structure\n", norm);

}

//------------------------------------------------------------------------------
// note: the heat power vector is *** NOT *** assembled

void StructExc::sendHeatPower(DistVec<double>& P) 
{

  double norm = 0.0;

  if (numStrNodes) {
    double (*buf)[1] = reinterpret_cast<double (*)[1]>(buffer);
#pragma omp parallel for reduction (+: norm)
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      Vec<double>& p = P(iSub);
      norm += p*p * pscale*pscale;
      matchNodes[iSub]->send(pscale, reinterpret_cast<double (*)[1]>(p.data()), buf);
    }

    for (int iCpu=0; iCpu<numStrCPU; ++iCpu) {
      if (numStrNodes[iCpu][0] > 0) {
	int size = numStrNodes[iCpu][0];
	double* localBuffer = buffer + numStrNodes[iCpu][1];
	strCom->sendTo(iCpu, HEATPOWER_TAG, localBuffer, size);
      }
    }
    strCom->waitForAllReq();
  }

  com->globalSum(1, &norm);

  com->printf(7, "Sent fluid heat power=%e to the structure\n", norm);

}

//------------------------------------------------------------------------------
void StructExc::getMdFreq(int &nf, double *&f)
{

  double buffer[2000];
  if (strCom->cpuNum() == 0)
    strCom->recFrom(1200, buffer, 1000);

  com->broadcast(1000, buffer);

  nf = int(buffer[0]);
  f = new double[nf];

  for(int i=0; i<nf; ++i)
    f[i] = buffer[i+1];
}

//------------------------------------------------------------------------------

void StructExc::getMdStrDisp(int id, DistSVec<double,3> &X0,
                             DistSVec<double,3> &X, DistSVec<double,3> &dX)
{
  double norms[2] = {0.0, 0.0};

  DistSVec<double,3> Xdot( X.info() );

  dX = 0.0;
  Xdot = 0.0;

  if (numStrNodes) {

    for (int iCpu=0; iCpu<numStrCPU; ++iCpu) {
      if (numStrNodes[iCpu][0] > 0) {
        int size = bufsize * numStrNodes[iCpu][0];
        double *localBuffer = buffer + bufsize * numStrNodes[iCpu][1];
        strCom->recFrom(iCpu, 1201 + id, localBuffer, size);
      }
    }

    double (*disp)[2][3] = reinterpret_cast<double (*)[2][3]>(buffer);

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {

      double (*x0)[3] = X0.subData(iSub);
      double (*x)[3] = X.subData(iSub);
      double (*xdot)[3] = Xdot.subData(iSub);
      double (*dx)[3] = dX.subData(iSub);

      double locNorms[2];

      matchNodes[iSub]->getDisplacement(algNum, dt, oolscale, oovscale, X.getMasterFlag(iSub),
                                        disp, x0, x, xdot, dx, locNorms);

#pragma omp critical
      norms[0] += locNorms[0];
#pragma omp critical
      norms[1] += locNorms[1];

    }

  }

  com->globalSum(2, norms);

  com->fprintf(stdout, "Received disp=%e and vel=%e from the structure \n", norms[0], norms[1]);
}
