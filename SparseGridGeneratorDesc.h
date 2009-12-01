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
  LocalRiemannGfmparGasJWL *lriemannGasJwl;

  void memberFunctionTest(double *in, double *res, double *parameters){ 
    res[0] = in[0]*in[0]*in[1]*in[1];
  }

public:

  SparseGridGeneratorDesc(IoData &ioData){
    refVal = new RefVal(ioData.ref.rv);
    varFcn = createVarFcn(ioData);
    lriemann = new LocalRiemann(varFcn);
    lriemannGasJwl = new LocalRiemannGfmparGasJWL(varFcn, ioData);

    if(varFcn->getType() == VarFcn::GAS) fprintf(stdout, "GAS VarFcn\n");
    else if(varFcn->getType() == VarFcn::LIQUID) fprintf(stdout, "LIQUID VarFcn\n");
    else if(varFcn->getType() == VarFcn::JWL) fprintf(stdout, "JWL VarFcn\n");
    else if(varFcn->getType() == VarFcn::GASINGAS) fprintf(stdout, "GASINGAS VarFcn\n");
    else if(varFcn->getType() == VarFcn::GASINLIQUID) fprintf(stdout, "GASINLIQUID VarFcn\n");
    else if(varFcn->getType() == VarFcn::LIQUIDINLIQUID) fprintf(stdout, "LIQUIDINLIQUID VarFcn\n");
    else if(varFcn->getType() == VarFcn::JWLINGAS) fprintf(stdout, "JWLINGAS VarFcn\n");
    else if(varFcn->getType() == VarFcn::JWLINJWL) fprintf(stdout, "JWLINJWL VarFcn\n");
    else fprintf(stdout, "no Type for VarFcn\n");
  }
  ~SparseGridGeneratorDesc(){
    delete refVal;
    delete varFcn;
    delete lriemann;
    //delete lriemannGasJwl;
  }

  void tabulate(IoData & ioData){
    srand(time(NULL));

    if(varFcn->getType() == VarFcn::JWL){

      double *refIn = new double[2]; double *refOut = new double[1];
      refIn[0] = ioData.ref.rv.density;
      refIn[1] = ioData.ref.rv.entropy;
      refOut[0] = ioData.ref.rv.velocity;
      fprintf(stdout, "SparseGridGeneratorDesc::tabulate::adim SparseGrid %e %e %e\n", refIn[0], refIn[1], refOut[0]);

      lriemann->setReferenceDensity(ioData.eqs.fluidModel.jwlModel.rhoref);
      fprintf(stdout, "### SparseGridGeneratorDesc::tabulate -- 1\n");
      double *parameters = new double;
      parameters[0] = 1.0;
      if(true){
        SparseGrid sparseGrid(ioData.mf.sparseGrid, parameters);
        sparseGrid.scaleGrid(refIn, refOut); // division operation
        fprintf(stdout, "### SparseGridGeneratorDesc::tabulate -- 2\n");
        //sparseGrid.tabulate(functionTest);
        //sparseGrid.tabulate(&SparseGridGeneratorDesc::memberFunctionTest,*this);
        sparseGrid.tabulate(&LocalRiemann::riemannInvariantGeneral2ndOrder,*lriemann);
        sparseGrid.printToFile(refIn, refOut); // multiplication operation
      }
      fprintf(stdout, "### SparseGridGeneratorDesc::tabulate -- 3\n");

      SparseGrid sparseGridCopy;
      sparseGridCopy.readFromFile();
      sparseGridCopy.scaleGrid(refIn, refOut);
      int numTest = 10;
      double **output = new double *[numTest+1];
      double **coord = new double *[numTest+1];
      fprintf(stdout, "range1 = %e %e\n", ioData.mf.sparseGrid.range[0][0],ioData.mf.sparseGrid.range[0][1]);
      fprintf(stdout, "range2 = %e %e\n", ioData.mf.sparseGrid.range[1][0],ioData.mf.sparseGrid.range[1][1]);
      for(int iTest=0; iTest<numTest+1; iTest++){
        output[iTest] = new double[1];
        coord[iTest]  = new double[2];
      }
      for(int iTest=0; iTest<numTest; iTest++){
        coord[iTest][0] = (ioData.mf.sparseGrid.range[0][0] + static_cast<double>(rand())/RAND_MAX * (ioData.mf.sparseGrid.range[0][1] - ioData.mf.sparseGrid.range[0][0]))/refIn[0];
        coord[iTest][1] = (ioData.mf.sparseGrid.range[1][0] + static_cast<double>(rand())/RAND_MAX * (ioData.mf.sparseGrid.range[1][1] - ioData.mf.sparseGrid.range[1][0]))/refIn[1];
        fprintf(stdout, "coord[%d] = (%e %e)\n", iTest, coord[iTest][0], coord[iTest][1]);
      }
      coord[10][0] = 1.007694e-01;
      coord[10][1] = 5.890818e+00;
      sparseGridCopy.interpolate(numTest+1, coord, output);
      double *exact = new double;
      for(int iTest=0; iTest<numTest+1; iTest++){
        //functionTest(coord[iTest],exact,parameters);
        //memberFunctionTest(coord[iTest],exact,parameters);
        lriemann->riemannInvariantGeneral2ndOrder(coord[iTest],exact,parameters);
        fprintf(stdout, "interpolation/exact output is %e/%e and relative error is %e\n", 
                output[iTest][0], exact[0], (exact[0]-output[iTest][0])/exact[0]);
      }

      for(int iTest=0; iTest<numTest+1; iTest++){
        delete [] output[iTest];
        delete [] coord[iTest];
      }
      delete exact; delete [] refIn; delete [] refOut; delete parameters; delete [] output; delete [] coord;

    }else if(varFcn->getType() == VarFcn::JWLINGAS){
      lriemannGasJwl->setReferenceDensity(ioData.eqs.fluidModel2.jwlModel.rhoref);
      double *parameters = new double;
      parameters[0] = 1.0;
      SparseGrid sparseGrid(ioData.mf.sparseGrid, parameters);
      //double *refIn = new double[5]; double *refOut = new double[1]; // 1output
      double *refIn = new double[5]; double *refOut = new double[2]; // 2outputs
      refIn[0] = ioData.ref.rv.density;
      refIn[1] = ioData.ref.rv.pressure;
      refIn[2] = ioData.ref.rv.density;
      refIn[3] = ioData.ref.rv.pressure;
      refIn[4] = ioData.ref.rv.velocity;
      refOut[0] = ioData.ref.rv.density;
      refOut[1] = ioData.ref.rv.density; // 2outputs
      fprintf(stdout, "refIn are %e %e %e\n", refIn[0],refIn[4],refIn[1]);


      /*double *testcoord = new double[6];
      double *testres   = new double[1];
      testcoord[0] = 1.039634e-01;
      testcoord[1] = 4.213924e-02;
      testcoord[2] = 1.073171e+00;
      testcoord[3] = 2.109069e+01;
      testcoord[4] = 1.314423e+01;
      lriemannGasJwl->eriemanngj_wrapper(testcoord,testres,parameters);
      exit(1);*/



      if(true){
        sparseGrid.scaleGrid(refIn, refOut); // division operation
        sparseGrid.tabulate(&LocalRiemannGfmparGasJWL::eriemanngj_wrapper,*lriemannGasJwl);
        sparseGrid.printToFile(refIn, refOut); // multiplication operation
      }

      SparseGrid sparseGridCopy;
      sparseGridCopy.readFromFile();
      sparseGridCopy.scaleGrid(refIn, refOut);
      // Testing the grid randomly
      int numTest = 50;
      double **output = new double *[numTest];
      double **coord = new double *[numTest];
      for(int iTest=0; iTest<numTest; iTest++){
        //output[iTest] = new double[1]; // 1output
        output[iTest] = new double[2]; // 2outputs
        coord[iTest]  = new double[5];
      }
      for(int iTest=0; iTest<numTest; iTest++){
        for(int iInput=0; iInput<5; iInput++)
          coord[iTest][iInput] = (ioData.mf.sparseGrid.range[iInput][0] + static_cast<double>(rand())/RAND_MAX * (ioData.mf.sparseGrid.range[iInput][1] - ioData.mf.sparseGrid.range[iInput][0]))/refIn[iInput];
        fprintf(stdout, "coord[%d] = (%e %e %e %e %e)\n", iTest, coord[iTest][0], coord[iTest][1], coord[iTest][2], coord[iTest][3], coord[iTest][4]);
      }
      sparseGridCopy.interpolate(numTest, coord, output);
      double *exact = new double;
      for(int iTest=0; iTest<numTest; iTest++){
        lriemannGasJwl->eriemanngj_wrapper(coord[iTest],exact,parameters);
        //fprintf(stdout, "interpolation/exact output is %e/%e and relative error is %e\n", 
        //        output[iTest][0], exact[0], (exact[0]-output[iTest][0])/exact[0]);  // 1output
        fprintf(stdout, "interpolation/exact output is (%e/%e, %e/%e) and relative error is (%e, %e)\n", 
                output[iTest][0], exact[0], output[iTest][1], exact[1], (exact[0]-output[iTest][0])/exact[0], (exact[1]-output[iTest][1])/exact[1]); // 2outputs
      }
      /*coord[0][0] = 6.6e-1;
      coord[0][1] = 2.65e-2;
      coord[0][2] = 3.2456;
      coord[0][3] = 2.4234e-1;
      coord[0][4] = 6.043135;
      coord[0][5] = 3.56681e1;
      fprintf(stdout, "coord[perso] = (%e %e %e %e %e %e)\n", coord[0][0], coord[0][1], coord[0][2], coord[0][3], coord[0][4], coord[0][5]);
      sparseGridCopy.interpolate(1, coord, output);
      lriemannGasJwl->eriemanngj_wrapper(coord[0],exact,parameters);
      fprintf(stdout, "interpolation/exact output is %e/%e and relative error is %e\n",output[0][0], exact[0], (exact[0]-output[0][0])/exact[0]);*/

      for(int iTest=0; iTest<numTest; iTest++){
        delete [] output[iTest];
        delete [] coord[iTest];
      }
      delete exact; delete [] output; delete [] coord;
      //delete parameters; delete [] refIn; delete [] refOut;

      if(true){
      // Testing the grid linearly...
      fprintf(stdout, "\n\n... Testing the grid linearly ...\n");
      numTest = 3125; //5 points in each 5 direction
      output = new double *[numTest];
      coord = new double *[numTest];
      for(int iTest=0; iTest<numTest; iTest++){
        //output[iTest] = new double[1]; // 1output
        output[iTest] = new double[2]; // 2outputs
        coord[iTest]  = new double[5];
      }
      int nIntervals = 5;
      double spacing[5]; double dimCoord[5][nIntervals];
      for (int iInput=0; iInput<5; iInput++){
        spacing[iInput] = (ioData.mf.sparseGrid.range[iInput][1] - ioData.mf.sparseGrid.range[iInput][0]) / refIn[iInput] /nIntervals;
        dimCoord[iInput][0] = ioData.mf.sparseGrid.range[iInput][0]/refIn[iInput] + 0.5*spacing[iInput];
        for (int iIntervals=1; iIntervals<nIntervals; iIntervals++)
          dimCoord[iInput][iIntervals] = dimCoord[iInput][iIntervals-1] + spacing[iInput];
      }
      for (int iInput=0; iInput<5; iInput++){
        fprintf(stdout, "dimBounds[%d] = [%e %e]\n", iInput, ioData.mf.sparseGrid.range[iInput][0]/ refIn[iInput], ioData.mf.sparseGrid.range[iInput][1]/ refIn[iInput]);
        fprintf(stdout, "dimCoord[%d] = %e %e %e %e %e\n", iInput, dimCoord[iInput][0],dimCoord[iInput][1],dimCoord[iInput][2],dimCoord[iInput][3],dimCoord[iInput][4]);
      }

      for(int idim1=0; idim1<5; idim1++)
        for(int idim2=0; idim2<5; idim2++)
          for(int idim3=0; idim3<5; idim3++)
            for(int idim4=0; idim4<5; idim4++)
              for(int idim5=0; idim5<5; idim5++){
                coord[idim1+5*(idim2+5*(idim3+5*(idim4+5*idim5)))][0] = dimCoord[0][idim1];
                coord[idim1+5*(idim2+5*(idim3+5*(idim4+5*idim5)))][1] = dimCoord[1][idim2];
                coord[idim1+5*(idim2+5*(idim3+5*(idim4+5*idim5)))][2] = dimCoord[2][idim3];
                coord[idim1+5*(idim2+5*(idim3+5*(idim4+5*idim5)))][3] = dimCoord[3][idim4];
                coord[idim1+5*(idim2+5*(idim3+5*(idim4+5*idim5)))][4] = dimCoord[4][idim5];
              }
      for(int iTest=0; iTest<numTest; iTest++)
        fprintf(stdout, "coord[%d] = (%e %e %e %e %e)\n", iTest, coord[iTest][0], coord[iTest][1], coord[iTest][2], coord[iTest][3], coord[iTest][4]);
      sparseGridCopy.interpolate(numTest, coord, output);
      double **exactt = new double *[numTest];
      for(int iTest=0; iTest<numTest; iTest++){
        //exactt[iTest] = new double[1]; // 1output
        exactt[iTest] = new double[2]; // 1output
        lriemannGasJwl->eriemanngj_wrapper(coord[iTest],exactt[iTest],parameters);
      }
      for(int iTest=0; iTest<numTest; iTest++)
        //fprintf(stdout, "%d %e %e %e %e %e %e %e %e\n", iTest, coord[iTest][0], coord[iTest][1], coord[iTest][2], coord[iTest][3], coord[iTest][4], output[iTest][0], exactt[iTest][0], (exactt[iTest][0]-output[iTest][0])/exactt[iTest][0]); // 1output
        fprintf(stdout, "%d %e %e %e %e %e (%e,%e) (%e,%e) (%e,%e)\n", iTest, coord[iTest][0], coord[iTest][1], coord[iTest][2], coord[iTest][3], coord[iTest][4], output[iTest][0], output[iTest][1], exactt[iTest][0], exactt[iTest][1], (exactt[iTest][0]-output[iTest][0])/exactt[iTest][0],(exactt[iTest][1]-output[iTest][1])/exactt[iTest][1]); // 2outputs
      
      for(int iTest=0; iTest<numTest; iTest++){
        delete [] output[iTest];
        delete [] coord[iTest];
        delete [] exactt[iTest];
      }
      delete [] exactt; delete [] output; delete [] coord;
      }
      if(false){
      // Testing the grid on grid points!
      fprintf(stdout, "\n\n... Testing the grid on the grid points ...\n");
      numTest = 4; 
      output = new double *[numTest];
      coord = new double *[numTest];
      for(int iTest=0; iTest<numTest; iTest++){
        output[iTest] = new double[1];
        coord[iTest]  = new double[5];
      }
      coord[0][0] = (ioData.mf.sparseGrid.range[0][0] + 0.5 * (ioData.mf.sparseGrid.range[0][1] - ioData.mf.sparseGrid.range[0][0]))/refIn[0];
      coord[0][1] = (ioData.mf.sparseGrid.range[1][0] + 0.5 * (ioData.mf.sparseGrid.range[1][1] - ioData.mf.sparseGrid.range[1][0]))/refIn[1];
      coord[0][2] = (ioData.mf.sparseGrid.range[2][0] + 0.5 * (ioData.mf.sparseGrid.range[2][1] - ioData.mf.sparseGrid.range[2][0]))/refIn[2];
      coord[0][3] = (ioData.mf.sparseGrid.range[3][0] + 0.5 * (ioData.mf.sparseGrid.range[3][1] - ioData.mf.sparseGrid.range[3][0]))/refIn[3];
      coord[0][4] = (ioData.mf.sparseGrid.range[4][0] + 0.5 * (ioData.mf.sparseGrid.range[4][1] - ioData.mf.sparseGrid.range[4][0]))/refIn[4];

      coord[1][0] = (ioData.mf.sparseGrid.range[0][0] + 0.25 * (ioData.mf.sparseGrid.range[0][1] - ioData.mf.sparseGrid.range[0][0]))/refIn[0];
      coord[1][1] = (ioData.mf.sparseGrid.range[1][0] + 0.25 * (ioData.mf.sparseGrid.range[1][1] - ioData.mf.sparseGrid.range[1][0]))/refIn[1];
      coord[1][2] = (ioData.mf.sparseGrid.range[2][0] + 0.25 * (ioData.mf.sparseGrid.range[2][1] - ioData.mf.sparseGrid.range[2][0]))/refIn[2];
      coord[1][3] = (ioData.mf.sparseGrid.range[3][0] + 0.25 * (ioData.mf.sparseGrid.range[3][1] - ioData.mf.sparseGrid.range[3][0]))/refIn[3];
      coord[1][4] = (ioData.mf.sparseGrid.range[4][0] + 0.25 * (ioData.mf.sparseGrid.range[4][1] - ioData.mf.sparseGrid.range[4][0]))/refIn[4];

      coord[2][0] = (ioData.mf.sparseGrid.range[0][0] + 0.00 * (ioData.mf.sparseGrid.range[0][1] - ioData.mf.sparseGrid.range[0][0]))/refIn[0];
      coord[2][1] = (ioData.mf.sparseGrid.range[1][0] + 0.50 * (ioData.mf.sparseGrid.range[1][1] - ioData.mf.sparseGrid.range[1][0]))/refIn[1];
      coord[2][2] = (ioData.mf.sparseGrid.range[2][0] + 0.75 * (ioData.mf.sparseGrid.range[2][1] - ioData.mf.sparseGrid.range[2][0]))/refIn[2];
      coord[2][3] = (ioData.mf.sparseGrid.range[3][0] + 0.50 * (ioData.mf.sparseGrid.range[3][1] - ioData.mf.sparseGrid.range[3][0]))/refIn[3];
      coord[2][4] = (ioData.mf.sparseGrid.range[4][0] + 0.75 * (ioData.mf.sparseGrid.range[4][1] - ioData.mf.sparseGrid.range[4][0]))/refIn[4];

      coord[3][0] = (ioData.mf.sparseGrid.range[0][0] + 0.00 * (ioData.mf.sparseGrid.range[0][1] - ioData.mf.sparseGrid.range[0][0]))/refIn[0];
      coord[3][1] = (ioData.mf.sparseGrid.range[1][0] + 0.50 * (ioData.mf.sparseGrid.range[1][1] - ioData.mf.sparseGrid.range[1][0]))/refIn[1];
      coord[3][2] = (ioData.mf.sparseGrid.range[2][0] + 0.25 * (ioData.mf.sparseGrid.range[2][1] - ioData.mf.sparseGrid.range[2][0]))/refIn[2];
      coord[3][3] = (ioData.mf.sparseGrid.range[3][0] + 0.50 * (ioData.mf.sparseGrid.range[3][1] - ioData.mf.sparseGrid.range[3][0]))/refIn[3];
      coord[3][4] = (ioData.mf.sparseGrid.range[4][0] + 0.00 * (ioData.mf.sparseGrid.range[4][1] - ioData.mf.sparseGrid.range[4][0]))/refIn[4];

      sparseGridCopy.interpolate(numTest, coord, output);
      double **exacttt = new double *[numTest];
      for(int iTest=0; iTest<numTest; iTest++){
        exacttt[iTest] = new double[1];
        lriemannGasJwl->eriemanngj_wrapper(coord[iTest],exacttt[iTest],parameters);
      }
      for(int iTest=0; iTest<numTest; iTest++)
        fprintf(stdout, "%d %e %e %e %e %e %e %e %e\n", iTest, coord[iTest][0], coord[iTest][1], coord[iTest][2], coord[iTest][3], coord[iTest][4], output[iTest][0], exacttt[iTest][0], (exacttt[iTest][0]-output[iTest][0])/exacttt[iTest][0]);
      
      for(int iTest=0; iTest<numTest; iTest++){
        delete [] output[iTest];
        delete [] coord[iTest];
        delete [] exacttt[iTest];
      }
      delete [] exacttt; delete [] output; delete [] coord;
      }
      if(false){
      fprintf(stdout, "\n\n... Testing the grid points and the missing grid points ...\n");
      int coef[5] = {1,0,1,0,0};
      for(int testnumber=0; testnumber<13; testnumber++){
      fprintf(stdout, "\n\n... testnumber = %d\n", testnumber);
      if(testnumber==0){ coef[0]=2;coef[1]=0;coef[2]=0;coef[3]=0;coef[4]=0;}
      if(testnumber==1){ coef[0]=0;coef[1]=2;coef[2]=0;coef[3]=0;coef[4]=0;}
      if(testnumber==2){ coef[0]=0;coef[1]=0;coef[2]=2;coef[3]=0;coef[4]=0;}
      if(testnumber==3){ coef[0]=0;coef[1]=0;coef[2]=0;coef[3]=2;coef[4]=0;}
      if(testnumber==4){ coef[0]=0;coef[1]=0;coef[2]=0;coef[3]=0;coef[4]=2;}
      if(testnumber==5){ coef[0]=0;coef[1]=0;coef[2]=1;coef[3]=1;coef[4]=1;}
      if(testnumber==6){ coef[0]=0;coef[1]=1;coef[2]=1;coef[3]=0;coef[4]=0;}
      if(testnumber==7){ coef[0]=0;coef[1]=1;coef[2]=0;coef[3]=1;coef[4]=0;}
      if(testnumber==8){ coef[0]=0;coef[1]=1;coef[2]=0;coef[3]=0;coef[4]=1;}
      if(testnumber==9){ coef[0]=1;coef[1]=1;coef[2]=0;coef[3]=0;coef[4]=0;}
      if(testnumber==10){ coef[0]=1;coef[1]=0;coef[2]=1;coef[3]=0;coef[4]=0;}
      if(testnumber==11){ coef[0]=1;coef[1]=0;coef[2]=0;coef[3]=1;coef[4]=0;}
      if(testnumber==12){ coef[0]=1;coef[1]=0;coef[2]=0;coef[3]=0;coef[4]=1;}
      // Testing the grid points and the missing grid point
      int numDimPoints[5];
      for (int i=0; i<5; i++){
        if(coef[i]==0) numDimPoints[i] = 1;
        if(coef[i]>=1) numDimPoints[i] = pow(2,coef[i])+1;
      }
      int numTest = numDimPoints[0]*numDimPoints[1]*numDimPoints[2]*numDimPoints[3]*numDimPoints[4];
      output = new double *[numTest];
      coord = new double *[numTest];
      double **sgcoord = new double *[numTest];
      for(int iTest=0; iTest<numTest; iTest++){
        output[iTest] = new double[1];
        coord[iTest]  = new double[5];
        sgcoord[iTest]  = new double[5];
      }
      double *dimCoord[5];
      for (int i=0; i<5; i++){
        dimCoord[i] = new double[numDimPoints[i]];
        if(numDimPoints[i]==1) dimCoord[i][0] = 0.5;
        else{
          dimCoord[i][0] = 0.0;
          for (int iPts=1; iPts<numDimPoints[i]-1; iPts++)
            dimCoord[i][iPts] = dimCoord[i][iPts-1] + pow(2,-coef[i]);
          dimCoord[i][numDimPoints[i]-1] = 1.0;
        }
      }

      for (int i=0; i<5; i++){
        fprintf(stdout, "dimCoord[%d] = ", i);
        for (int iPts=0; iPts<numDimPoints[i]; iPts++) fprintf(stdout, "%e ", dimCoord[i][iPts]);
        fprintf(stdout, "\n");
      }

      for(int idim1=0; idim1<numDimPoints[0]; idim1++)
        for(int idim2=0; idim2<numDimPoints[1]; idim2++)
          for(int idim3=0; idim3<numDimPoints[2]; idim3++)
            for(int idim4=0; idim4<numDimPoints[3]; idim4++)
              for(int idim5=0; idim5<numDimPoints[4]; idim5++){
                int num = idim1 + numDimPoints[0]*(idim2 + numDimPoints[1]*(idim3 + numDimPoints[2]*(idim4 + numDimPoints[3]*idim5)));
                sgcoord[num][0] = dimCoord[0][idim1];
                sgcoord[num][1] = dimCoord[1][idim2];
                sgcoord[num][2] = dimCoord[2][idim3];
                sgcoord[num][3] = dimCoord[3][idim4];
                sgcoord[num][4] = dimCoord[4][idim5];

                coord[num][0] = (ioData.mf.sparseGrid.range[0][0] + sgcoord[num][0] * (ioData.mf.sparseGrid.range[0][1] - ioData.mf.sparseGrid.range[0][0]))/refIn[0];
                coord[num][1] = (ioData.mf.sparseGrid.range[1][0] + sgcoord[num][1] * (ioData.mf.sparseGrid.range[1][1] - ioData.mf.sparseGrid.range[1][0]))/refIn[1];
                coord[num][2] = (ioData.mf.sparseGrid.range[2][0] + sgcoord[num][2] * (ioData.mf.sparseGrid.range[2][1] - ioData.mf.sparseGrid.range[2][0]))/refIn[2];
                coord[num][3] = (ioData.mf.sparseGrid.range[3][0] + sgcoord[num][3] * (ioData.mf.sparseGrid.range[3][1] - ioData.mf.sparseGrid.range[3][0]))/refIn[3];
                coord[num][4] = (ioData.mf.sparseGrid.range[4][0] + sgcoord[num][4] * (ioData.mf.sparseGrid.range[4][1] - ioData.mf.sparseGrid.range[4][0]))/refIn[4];
              }
      for(int iTest=0; iTest<numTest; iTest++)
        fprintf(stdout, "coord[%d] = (%e %e %e %e %e)\n", iTest, coord[iTest][0], coord[iTest][1], coord[iTest][2], coord[iTest][3], coord[iTest][4]);
      sparseGrid.interpolate(numTest, coord, output);
      //sparseGridCopy.interpolate(numTest, coord, output);
      double **exactttt = new double *[numTest];
      for(int iTest=0; iTest<numTest; iTest++){
        exactttt[iTest] = new double[1];
        lriemannGasJwl->eriemanngj_wrapper(coord[iTest],exactttt[iTest],parameters);
      }
      for(int iTest=0; iTest<numTest; iTest++)
        fprintf(stdout, "%d %e %e %e %e %e %e %e %e\n", iTest, sgcoord[iTest][0], sgcoord[iTest][1], sgcoord[iTest][2], sgcoord[iTest][3], sgcoord[iTest][4], output[iTest][0], exactttt[iTest][0], (exactttt[iTest][0]-output[iTest][0])/exactttt[iTest][0]);
      
      for(int iTest=0; iTest<numTest; iTest++){
        delete [] output[iTest];
        delete [] coord[iTest];
        delete [] sgcoord[iTest];
        delete [] exactttt[iTest];
      }
      delete [] exactttt; delete [] output; delete [] coord; delete [] sgcoord;
      }
      }

      delete parameters; delete [] refIn; delete [] refOut;
    }else{
      fprintf(stdout, "### SparseGridGeneratorDesc::nothing done!\n");
    }
    fprintf(stdout, "### SparseGridGeneratorDesc::tabulate -- finished\n");
  }

//------------------------------------------------------------------------------

  VarFcn *createVarFcn(IoData &ioData){
    VarFcn *vf = 0;
    if(ioData.eqs.numPhase == 1){
    if(ioData.eqs.fluidModel.fluid == FluidModelData::GAS)
      vf = new VarFcnPerfectGasEuler3D(ioData);
    else if(ioData.eqs.fluidModel.fluid == FluidModelData::LIQUID)
      vf = new VarFcnWaterCompressibleEuler3D(ioData);
    else if(ioData.eqs.fluidModel.fluid == FluidModelData::JWL)
      vf = new VarFcnJWLEuler3D(ioData);
    }else
      vf = new VarFcnJWLInGasEuler3D(ioData);

    if(!vf){
      fprintf(stdout, "*** Error: no valid choice for the VarFcn\n");
      exit(1);
    }
    return vf;
  }
    

};

//------------------------------------------------------------------------------


#endif
