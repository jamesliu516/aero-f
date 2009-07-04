#ifndef _SPARSEGRIDGENERATOR_DESC_H_
#define _SPARSEGRIDGENERATOR_DESC_H_

#include <IoData.h>

class RefVal;
class VarFcn;

void functionTest(double *in, double *res){ 
  res[0] = in[0]*in[0]*in[1] + in[0] + 100*in[1] + in[0]*in[1]/2.0;
  //fprintf(stderr, "in functionTest -- res = %e\n", res);
}

//------------------------------------------------------------------------------

class SparseGridGeneratorDesc {

private:

  RefVal *refVal;
  VarFcn *varFcn;
  LocalRiemann *lriemann;


public:

  SparseGridGeneratorDesc(IoData &ioData){
    refVal = new RefVal(ioData.ref.rv);
    varFcn = createVarFcn(ioData);
    lriemann = new LocalRiemann();

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
    fprintf(stderr, "### SparseGridGeneratorDesc::tabulate -- 1\n");
    if(varFcn->getType() == VarFcn::JWL){
      fprintf(stderr, "### SparseGridGeneratorDesc::tabulate -- 2\n");
      SparseGrid<2,1> sparseGrid(ioData.mf.sparseGrid);
      fprintf(stderr, "### SparseGridGeneratorDesc::tabulate -- 3\n");
      //sparseGrid.tabulate(lriemann->tabulateRiemann);
      sparseGrid.tabulate(functionTest);
      sparseGrid.printToFile();
      fprintf(stderr, "### SparseGridGeneratorDesc::tabulate -- 4\n");

      SparseGrid<2,1> sparseGridCopy;
      sparseGridCopy.readFromFile();
      typedef double Output[1];
      typedef double Coord[2];
      Output *output = new Output[1];
      Coord *coord = new Coord[1];
      coord[0][0] = 7.8745;
      coord[0][1] = 8.11111;
      sparseGridCopy.interpolate(1, coord, output);
      fprintf(stdout, "interpolation output is %e\n", output[0][0]);
      double *res = new double;
      functionTest(coord[0],res);
      fprintf(stdout, "exact output is %e\n", res[0]);

    }else{
      fprintf(stderr, "### SparseGridGeneratorDesc::tabulate -- 5\n");
    }
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
