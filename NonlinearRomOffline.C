#include <cstdio>
#include <cmath>
#include <sys/time.h>
#include <algorithm>
#include <cstdlib>
#include <Domain.h>
#include <NonlinearRomOffline.h>
#include <IoData.h>
#include <NonlinearRom.h>
#include <DistVector.h>
#include <VectorSet.h>
#include <MatVecProd.h>
#include <Timer.h>
#include <DistBcData.h>
#include <DistGeoState.h>
#include <DistTimeState.h>
#include <PostOperator.h>
#include <ParallelRom.h>


template <int dim>
NonlinearRomOfflineSolver<dim>::NonlinearRomOfflineSolver(Communicator *_com, IoData &_ioData, Domain &dom) :
          domain(dom), Xref(dom.getNodeDistInfo()), controlVol(dom.getNodeDistInfo()) {

 com = _com;
 ioData = &_ioData;
 geoState = 0;

}

//-------------------------------------------------------------------------------

template <int dim>
void NonlinearRomOfflineSolver<dim>::solve()  {

 // set up Timers
 Timer *modalTimer = domain.getTimer();
 double t0;

 const char *snapsFile = tInput->snapFile;

 if (ioData->problem.alltype == ProblemData::_ROB_CONSTRUCTION_ && snapsFile){
   t0 = modalTimer->getTime();
   NonlinearRom<dim> NonlinearRom(com,*ioData,domain);
   NonlinearRom.createAndAnalyzeDatabase();
   modalTimer->addPodConstrTime(t0);
 }
 else if (ioData->problem.alltype == ProblemData::_NONLINEAR_ROM_PREPROCESSING_){
   geoState = new DistGeoState(*ioData, &domain);
   geoState->setup1(tInput->positions, &Xref, &controlVol);
   GnatPreprocessing<dim> gappy(com,*ioData,domain,geoState);
   gappy.buildReducedModel();
 }
 else if (ioData->problem.alltype == ProblemData::_NONLINEAR_ROM_PREPROCESSING_STEP_1_){
   geoState = new DistGeoState(*ioData, &domain);
   geoState->setup1(tInput->positions, &Xref, &controlVol);
   GnatPreprocessingStep1<dim> gappy(com,*ioData,domain,geoState);
   gappy.buildReducedModel();
 }
 else if (ioData->problem.alltype == ProblemData::_NONLINEAR_ROM_PREPROCESSING_STEP_2_){
   geoState = new DistGeoState(*ioData, &domain);
   geoState->setup1(tInput->positions, &Xref, &controlVol);
   GnatPreprocessingStep2<dim> gappy(com,*ioData,domain,geoState);
   gappy.buildReducedModel();
 }
 else if (ioData->problem.alltype == ProblemData::_SURFACE_MESH_CONSTRUCTION_){
   geoState = new DistGeoState(*ioData, &domain);
   geoState->setup1(tInput->positions, &Xref, &controlVol);
   SurfMeshGen<dim> surfMeshBuilder(com,*ioData,domain,geoState);
   surfMeshBuilder.buildReducedModel();
 }
 else if (ioData->problem.alltype == ProblemData::_SAMPLE_MESH_SHAPE_CHANGE_){
   //KTC: CHANGE!!!
   geoState = new DistGeoState(*ioData, &domain);
   geoState->setup1(tInput->positions, &Xref, &controlVol);
   ReducedMeshShapeChanger<dim> reducedMeshShapeChanger(com,*ioData,domain,geoState);
   reducedMeshShapeChanger.buildReducedModel();
 }

}

//---------------------------------------------------------------------------------------
