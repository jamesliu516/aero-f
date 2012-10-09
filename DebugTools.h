// DebugTools.h
// Class with various tools for debugging
#ifndef _DEBUG_TOOLS_H_
#define _DEBUG_TOOLS_H_

#define GCCBACKTRACE
#ifdef GCCBACKTRACE
#include <unistd.h>
#include <execinfo.h>
#include <cstdio>
#include <cstdlib>
#endif


class DebugTools {

  public:

  template <int N,class Oper>
  static void CheckJacobian(const double jac[N*N], const double u[N],
                     const double f[N],const Oper& op, const char* name, double eps = 1.0e-7) {

    fprintf(stderr,"Checking jacobian for %s (output as (calc,act))\n", name);
    
    double fp[N], up[N];
    double Jp[N*N];
   
    for (int j = 0; j < N; ++j) {

      for (int k = 0; k < N; ++k)
        up[k] = u[k]*(1.0 + ((j!=k)?0.0:eps)) + ((j==k)?1.0e-7:0.0);
 
      op.Compute(up,fp);
      for (int k = 0; k < N; ++k)
        Jp[k*N+j] = (fp[k]-f[k])/(up[j]-u[j]);
    }

    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j)
        fprintf(stderr,"[%e %e] ", jac[i*N+j], Jp[i*N+j]);
      fprintf(stderr,"\n");
    }
  }

  static void PrintBacktrace() {

#ifdef GCCBACKTRACE
     void *array[100];
     size_t size;
     char **strings;
     size_t i;
     
     size = backtrace (array, 100);
     strings = backtrace_symbols (array, size);
     
     fprintf (stderr,"Backtrace from current location:\n"
                     "-----------------------------\n"
                     "Obtained %zd stack frames.\n", size);
     
     for (i = 0; i < size; i++)
        fprintf (stderr,"%s\n", strings[i]);
     
     free (strings);
#endif
   }
  static bool TryWaitForDebug() {

#ifdef AEROF_MPI_DEBUG
    int my_pid;// = getpid();
    MPI_Comm_rank(MPI_COMM_WORLD, &my_pid);
    FILE* file = fopen("aerofdebug","r"); 
    char fn[256];
    sprintf(fn, ".aerofdebug.pid.%d",getpid());
    std::cout << "PID " << getpid() << " is MPI rank " << my_pid << std::endl;
    if (!file)
      return false; 
    bool debug_process = false;
    while (!feof(file)) {
      int p;
      fscanf(file, "%d",&p);
      if (p == my_pid) {
        debug_process = true;
        break;
      }
    }
 
    volatile int wait = (debug_process?1:0);
    if (wait) {
 
      std::cout << "Process ID for rank " << my_pid << " is " << getpid() << std::endl;
    }

    while (wait) {

      sleep(5);
    }

    return debug_process;    
#else
    return false;
#endif
  }

  static void SpitRank() {

#ifdef AEROF_MPI_DEBUG

    int rnk;
    MPI_Comm_rank(MPI_COMM_WORLD,&rnk);
    std::cout << "Rank is " << rnk << std::endl;
#endif
  }

  template <class Scalar,int dim>
  static void PrintElement(const char* tag, DistSVec<Scalar,dim>& V,int rank,int iSub, int i) {

    int rnk;
    MPI_Comm_rank(MPI_COMM_WORLD,&rnk);
    if (rnk == rank) {

      std::cout << tag << "(" << iSub << ", " << i << ") = {";
      for (int k = 0; k < dim; ++k) 
        std::cout << V(iSub)[i][k] << ((k < dim-1)?", ":"}");
      std::cout << std::endl;
    }
  }
};

#endif
