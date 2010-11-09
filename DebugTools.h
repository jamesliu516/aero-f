// DebugTools.h
// Class with various tools for debugging

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

};
