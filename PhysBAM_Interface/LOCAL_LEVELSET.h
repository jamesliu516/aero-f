#ifndef __LOCAL_LEVELSET__
#define __LOCAL_LEVELSET__
#include<stdio.h>

class LocalLevelSet {
public:
  //NOTE: trId starts from 0.
  virtual double getPhi(int trId, double xi1, double xi2, bool* hasCracked=0) {
    fprintf(stderr,"ERROR: function getPhi is not implemented in LocalLevelSet!\n");
    return 0; 
  }
}; 

#endif
