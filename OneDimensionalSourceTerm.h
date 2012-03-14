#ifndef _ONE_DIMENSIONAL_SOURCE_TERM_H_
#define _ONE_DIMENSIONAL_SOURCE_TERM_H_

#include <cmath>

#include <DenseMatrixOps.h>

#include <sstream>

class OneDimensionalSourceTerm {

  struct MyLU {
    
    double* a;
    int* index;
  };

 public:
 
  OneDimensionalSourceTerm() { }

  ~OneDimensionalSourceTerm() { }

  int factorial(int i) { 
    return (i == 0 ? 1 : i*factorial(i-1)); 
  }

  void initialize(double _alpha, int _order, SVec<double,1>& X) {
    
    alpha = _alpha;
    order = _order;
    ludata = new MyLU[X.size()];
    int j;
    for (int i = 0; i < X.size(); ++i) {
      
      j = i-order/2;
      if (j < 0)
	j = 0;
      else if (j+order > X.size())
	j = X.size()-order;

      MyLU lu = {new double[order*order], new int[order]};
      for (int k = 0; k < order; ++k) {
	for (int l = 0; l < order; ++l) {
	  if (l == 0)
	    lu.a[k*order+l] = 1.0 / factorial(l);
	  else
	    lu.a[k*order+l] = pow( X[j+k][0]-X[i][0] , l) / factorial(l);
	}
      }
      DenseMatrixOp<double,5,5>::ludec(lu.a, lu.index,1.0, order);
      ludata[i] = lu;
    }
  }

  int binomialTerm(int n, int k) {
    return factorial(n)/(factorial(n-k)*factorial(k));
  }

  double computeIntegralTerm(double r1, double r0,double ri, int n) {
    
    double I = 0.0;
    if (r0 > 1.0e-8) {
      //std::cout << "ri = " << ri << " r1 = " << r1 <<  std::endl;
      I = pow(-ri, n)*log(r1/r0);
    }
    for (int k = 1; k <= n; ++k) {
      I += 1.0/k*binomialTerm(n,k)*(k==n?1.0 : pow(-ri,n-k))*(pow(r1,k)-pow(r0,k));
    }
    return I*alpha;
  }

  template <class FluxF,int dim>
    void compute(FluxF& f, SVec<double,dim>& V, SVec<double,dim>& F, SVec<double,1>& X, SVec<double,1>& Y, Vec<int>& fluidId) {
    
    double* local = new double[dim*V.size()];
    double* derivs = new double[order];

    int j;
    for (int i = 0; i < V.size(); ++i) {
      f.compute(V[i], local+i*dim,i);
    }

    //std::stringstream dummy;

    for (int i = 0; i < V.size(); ++i) {
      
      j = i-order/2;
      int fid = fluidId[i];
      if (j < 0)
	j = 0;
      else if (j+order > X.size())
	j = X.size()-order;

      MyLU& lu = ludata[i];
      for (int k = 0; k < dim; ++k) {
	for (int l = 0; l < order; ++l) {
          if (fluidId[j+l] == fid)
	    derivs[l] = local[(j+l)*dim+k];
          else
            derivs[l] = local[i*dim+k];
	}
	
	DenseMatrixOp<double,5,5>::ludfdbksb(lu.a, lu.index, derivs, order);
	//std::cout << derivs[0] << " " << factorial(0) << std::endl;
	
	double term = 0.0;
	for (int l = 0; l < order; ++l) {
	  //std::cout << "hello" << std::endl;
	  F[i][k] += computeIntegralTerm(Y[i+1][0],Y[i][0],X[i][0],l)*derivs[l]/factorial(l);
	  //l = l;
	  //std::cout << order << std::endl;// exit(0);
	  //std::cout << term << " ";
	  //dummy << "hello" << Y[i+1][0] << " " << Y[i][0] << " " << X[i][0] << " " << derivs[l] << std::endl;
	}
        //std::cout << term << " ";
	//exit(0);
	F[i][k] += term;
	//std::cout << term << " ";
      }
    } 

    delete [] local;
    delete [] derivs;
  }
  
 private:

  double alpha; // 1 for cylindrical, 2 for spherical

  int order;

  MyLU* ludata;

};

#endif
