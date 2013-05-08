#include <TsInput.h>
#include <sys/time.h>
#include <cstdlib>
#include <math.h>

using std::stable_sort;

template<class VecType>  
KspBinaryOutput<VecType>::KspBinaryOutput(Communicator *_com, IoData *_ioData, Domain *_domain)  : 
com(_com), ioData(_ioData), domain(_domain)
{ 
  krylovFreqTime = ioData->output.rom.krylovOutputFreqTime;
  krylovFreqNewton = ioData->output.rom.krylovOutputFreqNewton;
  krylovEnergy = ioData->output.rom.krylovVectorEnergy;
  timeIt = domain->getTimeIt();
  newtonIt = domain->getNewtonIt();

  fileName = new char[strlen(ioData->output.rom.prefix) + strlen(ioData->output.rom.krylovVector) + 1];
  sprintf(fileName, "%s%s", ioData->output.rom.prefix, ioData->output.rom.krylovVector);
}

//----------------------------------------------------------------------------------

template<class VecType> 
KspBinaryOutput<VecType>::~KspBinaryOutput() 
{

  delete [] fileName;
 
}

//----------------------------------------------------------------------------------

// this struct is used in writeKrylovVectors
struct kspSortStruct {
  int kspIndex; // index of corresponding krylov Vector
  double energy; // distance to second closest cluster

  bool operator<(const kspSortStruct& a) const {
    return energy < a.energy;
  }
};

//----------------------------------------------------------------------------------

template<class VecType>
template<int dim>
void KspBinaryOutput<VecType>::writeKrylovVectors(VecSet<DistSVec<double, dim> >& kspVecs, Vec<double> kspCoords, int numVecs) {

  if ((*(timeIt)%krylovFreqTime==0) && (*(newtonIt)%krylovFreqNewton==0)) {
    Vec<double> kspCoordsTrunc(numVecs);
    for (int iVec=0; iVec<numVecs; ++iVec) kspCoordsTrunc[iVec] = kspCoords[iVec];
    kspCoordsTrunc *= (1/kspCoordsTrunc.norm());
    kspSortStruct* kspIndexedCoords = new kspSortStruct[numVecs]; 
    for (int iVec=0; iVec<numVecs; ++iVec) {
      kspIndexedCoords[iVec].kspIndex = iVec;
      kspIndexedCoords[iVec].energy = pow(kspCoordsTrunc[iVec],2);
    }
    sort(kspIndexedCoords, kspIndexedCoords+numVecs); //ascending order
    double cumEnergy = 0;
    int vecsOutput = 0;
    while ((cumEnergy<krylovEnergy)&&(vecsOutput<numVecs)) {
      domain->writeVectorToFile(fileName, *(domain->getKrylovStep()), *(domain->getNewtonTag()), kspVecs[kspIndexedCoords[numVecs-vecsOutput-1].kspIndex]);
      cumEnergy += kspIndexedCoords[numVecs-vecsOutput-1].energy; 
      ++(*(domain->getKrylovStep()));
      ++vecsOutput;
      com->fprintf(stdout, "vecsOutput %d, energy %e, cumEnergy %e\n", vecsOutput, kspIndexedCoords[numVecs-vecsOutput-1].energy, cumEnergy);
    }
    delete [] kspIndexedCoords;
  }
}

//----------------------------------------------------------------------------------

template<class VecType>
template<int dim>
void KspBinaryOutput<VecType>::writeKrylovVectors(VecSet<DistSVec<bcomp, dim> >& kspVecs, Vec<bcomp> kspCoords, int numVecs) {
// do nothing
}

//----------------------------------------------------------------------------------

template<class VecType>
template<int dim, class Scalar>
void KspBinaryOutput<VecType>::writeKrylovVectors(VecSet<DistEmbeddedVec<Scalar, dim> >& kspVecs, Vec<Scalar> kspCoords, int numVecs) {
// do nothing
}

