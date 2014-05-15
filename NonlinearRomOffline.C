#include <cstdio>
#include <cmath>
#include <sys/time.h>
#include <algorithm>
#include <cstdlib>
#include <Domain.h>
#include <NonlinearRomOffline.h>
#include <IoData.h>
#include <NonlinearRom.h>
#include <NonlinearRomDatabaseConstruction.h>
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
NonlinearRomOfflineSolver<dim>::NonlinearRomOfflineSolver(Communicator *_com, IoData &_ioData, Domain &dom, GeoSource &_geoSource) :
          domain(dom), Xref(dom.getNodeDistInfo()), controlVol(dom.getNodeDistInfo()), geoSource(_geoSource) {

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


 if (ioData->problem.alltype == ProblemData::_NONLINEAR_ROM_OFFLINE_) {
   // create file system, cluster snapshots, construct ROBs
   if (ioData->romDatabase.nClusters > 0) {
     t0 = modalTimer->getTime();
     NonlinearRomDatabaseConstruction<dim> rom(com,*ioData,domain,geoSource);
     rom.constructDatabase();
     modalTimer->addPodConstrTime(t0);
   }
   // perform GNAT preprocessing (probably for snapshot collection method 0)
   const char *gnatPrefix = ioData->romDatabase.files.gnatPrefix;
   const char *sampledStateBasisName = ioData->romDatabase.files.sampledStateBasisName;
   if (strcmp(gnatPrefix,"")!=0 || strcmp(sampledStateBasisName,"")!=0) {
     geoState = new DistGeoState(*ioData, &domain);
     geoState->setup1(ioData->input.positions, &Xref, &controlVol);
     GnatPreprocessing<dim> gappy(com,*ioData,domain,geoState);
     gappy.buildReducedModel();
   }
   const char *surfacePrefix = ioData->romDatabase.files.surfacePrefix;
   const char *surfaceStateBasisName = ioData->romDatabase.files.surfaceStateBasisName;
   if (strcmp(surfacePrefix,"")!=0 || strcmp(surfaceStateBasisName,"")!=0) {
     geoState = new DistGeoState(*ioData, &domain);
     geoState->setup1(ioData->input.positions, &Xref, &controlVol);
     SurfMeshGen<dim> surfMeshGen(com,*ioData,domain,geoState);
     surfMeshGen.buildReducedModel();
   }
 }
 else if (ioData->problem.alltype == ProblemData::_NONLINEAR_ROM_PREPROCESSING_){
   geoState = new DistGeoState(*ioData, &domain);
   geoState->setup1(ioData->input.positions, &Xref, &controlVol);
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
   modalTimer->print(this->domain.getStrTimer());

 if (geoState) delete geoState;
}

//---------------------------------------------------------------------------------------
