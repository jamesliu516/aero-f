#ifndef _DIST_INFO_H_
#define _DIST_INFO_H_

#include <Communicator.h>

//------------------------------------------------------------------------------

struct DistInfo {

  int numLocThreads;
  int numLocSub;
  int totLen;
  int *subLen;
  int *subOffset;
  int *subLenReg;
  int *subOffsetReg;
  bool *masterFlag;
  double *invNdWeight; // inverse of the number of subs touching a node

  int numGlobSub;
  int *locSubToGlobSub;

  Communicator *com;

  DistInfo(int _numLocThreads, int _numLocSub, int _numGlobSub, 
	   int *_locSubToGlobSub, Communicator *_com)
  { 
    numLocThreads   = _numLocThreads;
    numLocSub       = _numLocSub;
    numGlobSub      = _numGlobSub;
    locSubToGlobSub = _locSubToGlobSub;
    com             = _com;
    subLen          = new int[numLocSub]; 
    subOffset       = new int[numLocSub];
    subLenReg       = new int[numLocThreads]; 
    subOffsetReg    = new int[numLocThreads];

    masterFlag     = 0;
    invNdWeight    = 0;
  }

  void setLen(int sub, int len) { subLen[sub] = len; }

  int subSize(int sub) const { return (subLen) ? subLen[sub] : 0; } //HB

  bool* getMasterFlag(int iSub) const
  {
    return (masterFlag) ? masterFlag+subOffset[iSub] : 0; 
  }

  double* getInvWeight(int iSub) const
  {
    return (invNdWeight) ? invNdWeight+subOffset[iSub] : 0;
  }

  void finalize(bool makeFlag) 
  {
    subOffset[0] = 0;

    for (int iSub = 1; iSub < numLocSub; ++iSub)
      subOffset[iSub] = subOffset[iSub-1] + subLen[iSub-1];

    totLen = subOffset[numLocSub-1] + subLen[numLocSub-1];

    if (makeFlag) {
       masterFlag = new bool[totLen];
       invNdWeight= new double[totLen];
    } else {
       masterFlag = 0;
       invNdWeight= 0;
    }

    numLocThreads = numLocSub;
    subLenReg     = subLen;
    subOffsetReg  = subOffset;    

    /* 
    int load = totLen / numLocThreads;
    int remainder = totLen % numLocThreads;

    for (int iSub=0; iSub<numLocThreads; ++iSub)
      subLenReg[iSub] = load;

    for (int iSub=0; iSub<remainder; ++iSub)
      subLenReg[iSub] += 1;

    subOffsetReg[0] = 0;

    for (int iSub = 1; iSub < numLocThreads; ++iSub)
      subOffsetReg[iSub] = subOffsetReg[iSub-1] + subLenReg[iSub-1];
    */
  }

  void print(char* mssg=0) //HB
  {
    com->sync();
    if(mssg) com->fprintf(stderr," DistInfo::print of %s\n",mssg);
    for(int iCPU=0; iCPU<com->size(); iCPU++) { 
      if(iCPU==com->cpuNum()) {
        fprintf(stderr," DistInfo::print on CPU %d:\n",com->cpuNum()); fflush(stderr);
        for(int iSub=0; iSub<numLocSub; ++iSub) {
          fprintf(stderr," * subd %3d, subOffset[%3d] = %6d, subLen[%3d] = %6d\n",
                locSubToGlobSub[iSub],iSub,subOffset[iSub],iSub,subLen[iSub]); fflush(stderr);
        }
      }
      com->sync();
    }
  }

  ~DistInfo() 
  {
    /* HB: according to Thuan, this causes troubles when running Linearized/Modal analysis
    if(subLen)       { delete [] subLen;       subLen      = 0; }
    if(subOffset)    { delete [] subOffset;    subOffset   = 0; }
    if(subLenReg)    { delete [] subLenReg;    subLenReg   = 0; }
    if(subOffsetReg) { delete [] subOffsetReg; subOffsetReg= 0; }
    if(masterFlag)   { delete [] masterFlag;   masterFlag  = 0; }
    if(invNdWeight)  { delete [] invNdWeight;  invNdWeight = 0; }
    */
    // Adam 2010.07.29 
    // It is needed though to detect memory leaks. Otherwise Valgrind complains…
    // Adam 2010.08.04 (Update)
    // Someone has to spend some time on the Valgrind errors. The objects are not destroyed correctly at the end of the simulation. 
    // If the following lines are uncommented, the simulation get stuck forever into this destructor... Memory seems to be corrupted somehow and this is not a good sign.
    /*
    delete [] subLen; 
    delete [] invNdWeight;
    delete [] subOffset;
    delete [] subLenReg;
    delete [] subOffsetReg;
    delete [] masterFlag;
    */
    /*
    fprintf(stderr,"yep!1\n");
    fprintf(stderr,"%p %i \n",invNdWeight,totLen);
    fprintf(stderr,"yep!2\n");
    for(int i=0;i<totLen;++i) {fprintf(stderr,"%f ",invNdWeight[i]);}
    fprintf(stderr,"\nyep!3\n");
    */
  }

  // Copy Constructor and Assignement Operator are declared Private
  // Compiler will return an error if used.
private:
  DistInfo(const DistInfo &);
  DistInfo& operator=(const DistInfo &);

};

//------------------------------------------------------------------------------
#endif
