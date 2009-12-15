#ifndef _SPARSEGRIDGENERATOR_DESC_H_
#define _SPARSEGRIDGENERATOR_DESC_H_

#include <IoData.h>
#include <SparseGrid.h>
#include <SparseGridCluster.h>
#include <cmath>
#include <ctime>
#include <Communicator.h>

class RefVal;
class VarFcn;
class LocalRiemann;

//------------------------------------------------------------------------------
  void functionTest(double *in, double *res, double *parameters){ 
    res[0] = in[0]*in[0]*in[1]*in[1];
  }

//------------------------------------------------------------------------------
class SparseGridGeneratorDesc {

private:

  VarFcn *varFcn;
  LocalRiemannGfmparGasJWL *lriemannGasJwl;
  Communicator *com;

  void memberFunctionTest(double *in, double *res, double *parameters){ 
    res[0] = in[0]*in[0]*in[1]*in[1];
  }

public:

  SparseGridGeneratorDesc(IoData &ioData, Communicator *comm){
    com = comm;
    varFcn = createVarFcn(ioData);
    lriemannGasJwl = new LocalRiemannGfmparGasJWL(varFcn, NULL, MultiFluidData::RK2);

    if(varFcn->getType() == VarFcn::GAS) com->fprintf(stdout, "GAS VarFcn\n");
    else if(varFcn->getType() == VarFcn::LIQUID) com->fprintf(stdout, "LIQUID VarFcn\n");
    else if(varFcn->getType() == VarFcn::JWL) com->fprintf(stdout, "JWL VarFcn\n");
    else if(varFcn->getType() == VarFcn::GASINGAS) com->fprintf(stdout, "GASINGAS VarFcn\n");
    else if(varFcn->getType() == VarFcn::GASINLIQUID) com->fprintf(stdout, "GASINLIQUID VarFcn\n");
    else if(varFcn->getType() == VarFcn::LIQUIDINLIQUID) com->fprintf(stdout, "LIQUIDINLIQUID VarFcn\n");
    else if(varFcn->getType() == VarFcn::JWLINGAS) com->fprintf(stdout, "JWLINGAS VarFcn\n");
    else if(varFcn->getType() == VarFcn::JWLINJWL) com->fprintf(stdout, "JWLINJWL VarFcn\n");
    else com->fprintf(stdout, "no Type for VarFcn\n");
  }

  ~SparseGridGeneratorDesc(){
    delete varFcn;
    delete lriemannGasJwl;
    com = 0;
  }

//------------------------------------------------------------------------------

  void tabulate(IoData & ioData){
    srand(time(NULL));

    if(ioData.mf.riemannComputation == MultiFluidData::TABULATION2){

      double *refIn = new double[2]; double *refOut = new double[1];
      refIn[0] = ioData.ref.rv.density;
      refIn[1] = ioData.ref.rv.entropy;
      refOut[0] = ioData.ref.rv.velocity;
      com->fprintf(stdout, "refIn/refOut are %e %e / %e\n", refIn[0], refIn[1], refOut[0]);

      double *parameters = new double[2];
      parameters[0] = 1.0;
      parameters[1] = ioData.eqs.fluidModel.jwlModel.rhoref;

      if(false){
        SparseGrid sparseGrid(ioData.mf.sparseGrid, parameters, refIn, refOut);
        com->fprintf(stdout, "### SparseGridGeneratorDesc::tabulate -- 2\n");
        //sparseGrid.tabulate(functionTest);
        //sparseGrid.tabulate(&SparseGridGeneratorDesc::memberFunctionTest,*this);
        sparseGrid.tabulate(&LocalRiemannGfmparGasJWL::riemannInvariantGeneral2ndOrder_wrapper,*lriemannGasJwl);
        sparseGrid.printToFile(refIn, refOut, ioData.output.transient.sparseGrid);
      }
      com->fprintf(stdout, "### SparseGridGeneratorDesc::tabulate -- 3\n");

      SparseGrid sparseGridCopy;
      sparseGridCopy.readFromFile(refIn, refOut, ioData.output.transient.sparseGrid);
      int number = 5;
      sparseGridCopy.test(&LocalRiemannGfmparGasJWL::riemannInvariantGeneral2ndOrder_wrapper,*lriemannGasJwl, 2, &number, parameters);

      delete [] refIn; delete [] refOut; delete [] parameters;
    }else if(ioData.mf.riemannComputation == MultiFluidData::TABULATION5){
      double *parameters = NULL;
      double *refIn = new double[5]; double *refOut = new double[2]; // 2outputs
      refIn[0] = ioData.ref.rv.density;
      refIn[1] = ioData.ref.rv.pressure;
      refIn[2] = ioData.ref.rv.density;
      refIn[3] = ioData.ref.rv.pressure;
      refIn[4] = ioData.ref.rv.velocity;
      refOut[0] = ioData.ref.rv.density;
      refOut[1] = ioData.ref.rv.density; // 2outputs
      com->fprintf(stdout, "refIn are %e %e %e\n", refIn[0],refIn[4],refIn[1]);

/*
      SparseGrid sparseGrid(ioData.mf.sparseGrid, parameters, refIn, refOut);
      if(true){
        sparseGrid.tabulate(&LocalRiemannGfmparGasJWL::eriemanngj_wrapper,*lriemannGasJwl);
        sparseGrid.printToFile(refIn, refOut, ioData.output.transient.sparseGrid);
      }

      SparseGrid sparseGridCopy;
      sparseGridCopy.readFromFile(refIn, refOut, ioData.output.transient.sparseGrid);
      int numTest = 5;
      if(true) 
        sparseGridCopy.test(&LocalRiemannGfmparGasJWL::eriemanngj_wrapper,*lriemannGasJwl,1,&numTest, parameters);

      if(true)
        sparseGridCopy.test(&LocalRiemannGfmparGasJWL::eriemanngj_wrapper,*lriemannGasJwl,2,&numTest, parameters);
*/

      SparseGridCluster sgCluster;
      sgCluster.generate(ioData.mf.sparseGrid, parameters, &LocalRiemannGfmparGasJWL::eriemanngj_wrapper,*lriemannGasJwl, ioData.output.transient.sparseGrid, refIn, refOut, com);
      delete [] refIn; delete [] refOut; delete parameters;
    }else{
      com->fprintf(stdout, "### SparseGridGeneratorDesc::nothing done!\n");
    }
    com->fprintf(stdout, "### SparseGridGeneratorDesc::tabulate -- finished\n");
  }

//------------------------------------------------------------------------------

  VarFcn *createVarFcn(IoData &ioData){
    VarFcn *vf = 0;
    if(ioData.mf.riemannComputation == MultiFluidData::TABULATION2){
      vf = new VarFcnJWLEuler3D(ioData);
    }else if(ioData.mf.riemannComputation == MultiFluidData::TABULATION5){
      vf = new VarFcnJWLInGasEuler3D(ioData);
    }

    if(!vf){
      com->fprintf(stdout, "*** Error: no valid choice for the VarFcn\n");
      exit(1);
    }
    return vf;
  }
    

};

//------------------------------------------------------------------------------


#endif
