#ifndef _ERROR_HANDLER_H_
#define _ERROR_HANDLER_H_

#include <Communicator.h>
#include <vector>

struct ErrorHandler{

  enum Error {UNPHYSICAL = 0, SATURATED_LS = 1, BAD_RIEMANN = 2, REDUCE_TIMESTEP = 3, PRESSURE_CLIPPING = 4, REDO_TIMESTEP = 5, LARGE_VELOCITY = 6, RAPIDLY_CHANGING_PRESSURE = 7, SIZE = 8};
  enum Type {ALL=0, SOLVER = 1};
  int localErrors[SIZE];
  int globalErrors[SIZE];
  //int solverErrors[];
  std::vector<int> *solverErrors;
  Communicator *com;

  ErrorHandler(Communicator *comIn){
    com = comIn;
    int solverErrorsArray[] = {UNPHYSICAL,SATURATED_LS,BAD_RIEMANN,PRESSURE_CLIPPING,LARGE_VELOCITY, RAPIDLY_CHANGING_PRESSURE};
    solverErrors = new std::vector<int>(solverErrorsArray,solverErrorsArray + sizeof(solverErrorsArray)/sizeof(int));
    for (int i=0; i<SIZE; i++) localErrors[i]=0;
  }

  void reduceError(){ 
    for (int i=0; i<SIZE; i++) globalErrors[i]=localErrors[i];
    com->globalSum(SIZE,globalErrors);
    return;
  }

  void clearError(int type=ALL){
    if(type==ALL) for(int i=0; i<SIZE; i++) {globalErrors[i]=0; localErrors[i]=0; }
    if(type==SOLVER) for(int i=0; i<solverErrors->size(); i++) {globalErrors[solverErrors->at(i)]=0; localErrors[solverErrors->at(i)]=0;}
    return;
  }

  void printError(int type=ALL){
    char str[200];
    sprintf(str,"");
    if(type==ALL) for(int i=0; i<SIZE; i++) {sprintf(str,"%s%i, ",str,globalErrors[i]);}
    if(type==SOLVER) for(int i=0; i<solverErrors->size(); i++) {sprintf(str,"%s%i, ",str,globalErrors[solverErrors->at(i)]);}
    //com->printf(1,"%s\n",str);
    std::printf("%s\n",str); //Race condition?
    fflush(stdout);
    com->barrier();
    return;

  }

};


#endif
