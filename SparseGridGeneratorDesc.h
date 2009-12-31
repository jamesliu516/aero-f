#ifndef _SPARSEGRIDGENERATOR_DESC_H_
#define _SPARSEGRIDGENERATOR_DESC_H_

#include <IoData.h>
#include <SparseGrid.h>
#include <math.h>
#include <time.h>

class RefVal;
class VarFcn;
class LocalRiemann;

//------------------------------------------------------------------------------
  void functionTest(double *in, double *res, double *parameters){ 
    res[0] = in[0]*in[0]*in[1]*in[1];
  }

class SparseGridGeneratorDesc {

private:

  RefVal *refVal;
  VarFcn *varFcn;
  LocalRiemann *lriemann;

  void memberFunctionTest(double *in, double *res, double *parameters){ 
    res[0] = in[0]*in[0]*in[1]*in[1];
  }

public:

  SparseGridGeneratorDesc(IoData &ioData){
    refVal = new RefVal(ioData.ref.rv);
    varFcn = createVarFcn(ioData);
    lriemann = new LocalRiemann(varFcn);

    if(varFcn->getType() == VarFcn::GAS) fprintf(stderr, "GAS VarFcn\n");
    else if(varFcn->getType() == VarFcn::LIQUID) fprintf(stderr, "LIQUID VarFcn\n");
    else if(varFcn->getType() == VarFcn::JWL) fprintf(stderr, "JWL VarFcn\n");
    else if(varFcn->getType() == VarFcn::GASINGAS) fprintf(stderr, "GASINGAS VarFcn\n");
    else if(varFcn->getType() == VarFcn::GASINLIQUID) fprintf(stderr, "GASINLIQUID VarFcn\n");
    else if(varFcn->getType() == VarFcn::LIQUIDINLIQUID) fprintf(stderr, "LIQUIDINLIQUID VarFcn\n");
    else if(varFcn->getType() == VarFcn::JWLINGAS) fprintf(stderr, "JWLINGAS VarFcn\n");
    else if(varFcn->getType() == VarFcn::JWLINJWL) fprintf(stderr, "JWLINJWL VarFcn\n");
    else fprintf(stderr, "no Type for VarFcn\n");
  }
  ~SparseGridGeneratorDesc(){
    delete refVal;
    delete varFcn;
    delete lriemann;
  }

  void tabulate(IoData & ioData){
    srand(time(NULL));

    double *refIn = new double[2]; double *refOut = new double[1];
    refIn[0] = ioData.ref.rv.density;
    refIn[1] = ioData.ref.rv.entropy;
    refOut[0] = ioData.ref.rv.velocity;
    fprintf(stdout, "adim SparseGrid %e %e %e\n", refIn[0], refIn[1], refOut[0]);

    if(varFcn->getType() == VarFcn::JWL){

      lriemann->setReferenceDensity(ioData.eqs.fluidModel.jwlModel.rhoref);
      fprintf(stderr, "### SparseGridGeneratorDesc::tabulate -- 1\n");
      double *parameters = new double;
      parameters[0] = 1.0;
      SparseGrid sparseGrid(ioData.mf.sparseGrid, parameters);
      sparseGrid.scaleGrid(refIn, refOut); // division operation
      fprintf(stderr, "### SparseGridGeneratorDesc::tabulate -- 2\n");
      //sparseGrid.tabulate(functionTest);
      //sparseGrid.tabulate(&SparseGridGeneratorDesc::memberFunctionTest,*this);
      sparseGrid.tabulate(&LocalRiemann::riemannInvariantGeneral2ndOrder,*lriemann);
      sparseGrid.printToFile(refIn, refOut); // multiplication operation
      fprintf(stderr, "### SparseGridGeneratorDesc::tabulate -- 3\n");

      SparseGrid sparseGridCopy;
      sparseGridCopy.readFromFile();
      sparseGridCopy.scaleGrid(refIn, refOut);
      int numTest = 10;
      double **output = new double *[numTest];
      double **coord = new double *[numTest];
      fprintf(stdout, "range1 = %e %e\n", ioData.mf.sparseGrid.range[0][0],ioData.mf.sparseGrid.range[0][1]);
      fprintf(stdout, "range2 = %e %e\n", ioData.mf.sparseGrid.range[1][0],ioData.mf.sparseGrid.range[1][1]);
      for(int iTest=0; iTest<numTest; iTest++){
        output[iTest] = new double[1];
        coord[iTest]  = new double[2];
        coord[iTest][0] = (ioData.mf.sparseGrid.range[0][0] + static_cast<double>(rand())/RAND_MAX * (ioData.mf.sparseGrid.range[0][1] - ioData.mf.sparseGrid.range[0][0]))/refIn[0];
        coord[iTest][1] = (ioData.mf.sparseGrid.range[1][0] + static_cast<double>(rand())/RAND_MAX * (ioData.mf.sparseGrid.range[1][1] - ioData.mf.sparseGrid.range[1][0]))/refIn[1];
        fprintf(stdout, "coord[%d] = (%e %e)\n", iTest, coord[iTest][0], coord[iTest][1]);
      }
      sparseGridCopy.interpolate(numTest, coord, output);
      double *exact = new double;
      for(int iTest=0; iTest<numTest; iTest++){
        //functionTest(coord[iTest],exact,parameters);
        //memberFunctionTest(coord[iTest],exact,parameters);
        lriemann->riemannInvariantGeneral2ndOrder(coord[iTest],exact,parameters);
        fprintf(stdout, "interpolation/exact output is %e/%e and relative error is %e\n", 
                output[iTest][0], exact[0], (exact[0]-output[iTest][0])/exact[0]);
      }

    }else{
      fprintf(stderr, "### SparseGridGeneratorDesc::nothing done!\n");
    }
    fprintf(stderr, "### SparseGridGeneratorDesc::tabulate -- finished\n");
  }


  VarFcn *createVarFcn(IoData &ioData){
    VarFcn *vf = 0;
    if(ioData.eqs.fluidModel.fluid == FluidModelData::GAS)
      vf = new VarFcnPerfectGasEuler3D(ioData);
    else if(ioData.eqs.fluidModel.fluid == FluidModelData::LIQUID)
      vf = new VarFcnWaterCompressibleEuler3D(ioData);
    else if(ioData.eqs.fluidModel.fluid == FluidModelData::JWL)
      vf = new VarFcnJWLEuler3D(ioData);

    if(!vf){
      fprintf(stderr, "*** Error: no valid choice for the VarFcn\n");
      exit(1);
    }
    return vf;
  }
    

};

//------------------------------------------------------------------------------


#endif
