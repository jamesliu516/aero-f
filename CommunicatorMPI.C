#include <alloca.h>
#include <stdio.h>
#include <stdlib.h>

#include <Communicator.h>

#ifndef MPI_INTEGER
#define MPI_INTEGER MPI_INT
#endif

#include <complex>
using std::complex;

#ifdef USE_MPI
template<>
MPI_Datatype CommTrace<int>::MPIType = MPI_INTEGER;
template<>
MPI_Datatype CommTrace<float>::MPIType = MPI_FLOAT;
template<>
MPI_Datatype CommTrace<double>::MPIType = MPI_DOUBLE;
//template<>
//MPI_Datatype CommTrace<complex<double> >::MPIType = MPI_DOUBLE_COMPLEX;
template<>
MPI_Datatype CommTrace<complex<double> >::MPIType = MPI_DOUBLE;

template<>
int CommTrace<int>::multiplicity = 1;
template<>
int CommTrace<float>::multiplicity = 1;
template<>
int CommTrace<double>::multiplicity = 1;
template<>
int CommTrace<complex<double> >::multiplicity = 2;

static MPI_Request nullReq;
static MPI_Status nullStat;
#endif

//------------------------------------------------------------------------------

void initCommunication(int &argc, char **&argv)
{

#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif

}

//------------------------------------------------------------------------------

void closeCommunication()
{

#ifdef USE_MPI
  MPI_Finalize();
#endif

}

//------------------------------------------------------------------------------

/*
#ifdef USE_MPI
void exit(int status)
{

  MPI_Abort(MPI_COMM_WORLD, status);

}
#endif
*/

//------------------------------------------------------------------------------

Communicator::Communicator()
#ifdef USE_MPI
  : pendReq(nullReq), reqStatus(nullStat)
#endif
{

#ifdef USE_MPI
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &numCPU);
  MPI_Comm_rank(comm, &thisCPU);
  nPendReq = 0;
#else
  numCPU = 1;
  thisCPU = 0;
#endif

  timer = 0;
  maxverbose = 0;

}

//------------------------------------------------------------------------------

#ifdef USE_MPI
Communicator::Communicator(MPI_Comm c1)
  : pendReq(nullReq), reqStatus(nullStat)
{

  comm = c1;
  MPI_Comm_size(comm, &numCPU);
  MPI_Comm_rank(comm, &thisCPU);
  nPendReq = 0;

  timer = 0;
  maxverbose = 0;

}
#endif

//------------------------------------------------------------------------------
// note: color+1 is used in order to make the routine compatible with the 
// fortran communication library (which is restricted to 4 codes)

void Communicator::split(int color, int maxcolor, Communicator** c)
{

  int i;
  for (i=0; i<maxcolor; ++i)
    c[i] = 0;

#ifdef USE_MPI
  int rank;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm comm1;
  MPI_Comm_split(comm, color+1, rank, &comm1);
  c[color] = new Communicator(comm1);

  int* leaders = new int[maxcolor];
  int* newleaders = new int[maxcolor];
  for (i=0; i<maxcolor; ++i)
    leaders[i] = -1;

  int localRank;
  MPI_Comm_rank(comm1, &localRank);
  if (localRank == 0)
    leaders[color] = rank;

  MPI_Allreduce(leaders, newleaders, maxcolor, MPI_INTEGER, MPI_MAX, comm); // May be a prob

  for (i=0; i<maxcolor; ++i) {
    if (i != color && newleaders[i] >= 0) {
      int tag;
      if (color < i)
	tag = maxcolor*(color+1)+i+1;
      else
	tag = maxcolor*(i+1)+color+1;
      MPI_Comm comm2;
      MPI_Intercomm_create(comm1, 0, comm, newleaders[i], tag, &comm2);
      c[i] = new Communicator(comm2);
    }
  }

  if (leaders) 
    delete [] leaders;
  if (newleaders) 
    delete [] newleaders;
#else
  c[color] = this;
#endif

}

//------------------------------------------------------------------------------

int Communicator::remoteSize()
{

  int numRemote = 0;

#ifdef USE_MPI
  MPI_Comm_remote_size(comm, &numRemote);
#endif

  return numRemote;

}

//------------------------------------------------------------------------------

int Communicator::barrier()
{

  int ierr = 0;

#ifdef USE_MPI
  ierr = MPI_Barrier(comm); 
#endif

  return ierr;

}

//------------------------------------------------------------------------------


void Communicator::printf(int verbose, const char *format, ...)
{

  if (thisCPU == 0 && verbose <= maxverbose) {
    va_list args;
    va_start(args, format);
    vfprintf(stdout, format, args);
    ::fflush(stdout);
    va_end(args);
  }

}

//------------------------------------------------------------------------------

void Communicator::fprintf(FILE *file, const char *format, ...)
{

  if (thisCPU == 0) {
    va_list args;
    va_start(args, format);
    vfprintf(file, format, args);
    ::fflush(file);
    va_end(args);
  }

}

//------------------------------------------------------------------------------

void Communicator::waitForAllReq()
{

#ifdef USE_MPI
  // Just making sure that reqStatus has an appropriate length
  MPI_Status *safe = reqStatus+nPendReq;

  if (safe == 0)
    exit(1);

  int nSuccess = MPI_Waitall(nPendReq, pendReq+0, reqStatus+0);

  if (nSuccess) {
    fprintf(stderr, "*** Error: unexpected success number %d\n", nSuccess);
    exit(1);
  }

  nPendReq = 0;
#endif

}

//------------------------------------------------------------------------------
