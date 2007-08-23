#ifndef _TIME_STATE_H_
#define _TIME_STATE_H_

class TimeData;
class GeoState;

template<class Scalar> class Vec;
template<class Scalar, int dim> class SVec;
template<class Scalar, int dim> class GenMat;

//------------------------------------------------------------------------------

template<int dim>
class TimeState {

  TimeData &data;

  Vec<double> &dt;
  Vec<double> &idti;
  Vec<double> &idtv;
  SVec<double,dim> &Un;
  SVec<double,dim> &Unm1;
  SVec<double,dim> &Unm2;
  SVec<double,dim> &Rn;

public:

  TimeState(TimeData &, Vec<double> &, Vec<double> &, Vec<double> &, SVec<double,dim> &, 
	    SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &);
  ~TimeState() {}

  void add_dAW_dt(bool *, GeoState &, Vec<double> &, 
			  SVec<double,dim> &, SVec<double,dim> &);
  void add_dAW_dtLS(bool *, GeoState &, Vec<double> &, 
			  Vec<double> &, Vec<double> &, Vec<double> &, Vec<double> &);

  template<class Scalar, int neq>
  void addToJacobianNoPrec(bool *, Vec<double> &, GenMat<Scalar,neq> &, SVec<double,dim> &,
                     VarFcn *, int*);
  template<class Scalar, int neq>
  void addToJacobianNoPrecLocal(int, double, SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<class Scalar, int neq>
  void addToJacobianGasPrec(bool *, Vec<double> &, GenMat<Scalar,neq> &, SVec<double,dim> &,
                     VarFcn *, double, double, double, double, double, Vec<double> &, int*);
  template<class Scalar, int neq>
  void addToJacobianGasPrecLocal(int, double, double, double, double, double, double, double,
				 SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<class Scalar, int neq>
  void addToJacobianLiquidPrec(bool *, Vec<double> &, GenMat<Scalar,neq> &, SVec<double,dim> &,
                     VarFcn *, double, double, double, Vec<double> &, int*);
  template<class Scalar, int neq>
  void addToJacobianLiquidPrecLocal(int, double, VarFcn *, double, double, double, double,
				    SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<class Scalar, int neq>
  void addToH1(bool *, Vec<double> &, GenMat<Scalar,neq> &);

  template<class Scalar, int neq>
  void addToH1(bool *, Vec<double> &, GenMat<Scalar,neq> &, Scalar);

  template<class Scalar>
  void addToH2(bool *, VarFcn *, Vec<double> &, SVec<double,dim> &, GenMat<Scalar,dim> &);

  template<class Scalar>
  void addToH2(bool *, VarFcn *, Vec<double> &, SVec<double,dim> &,
               GenMat<Scalar,dim> &, Scalar);

  template<class Scalar>
  void addToH2LS(bool *, VarFcn *, Vec<double> &, SVec<double,dim> &, GenMat<Scalar,1> &);

  template<class Scalar>
  void addToH2Minus(bool *, VarFcn *, Vec<double> &, SVec<double,dim> &, GenMat<Scalar,dim> &);
                                                                                                                          
  void get_dW_dt(bool *, GeoState &, Vec<double> &, SVec<double,dim> &, SVec<double,dim> &);
                                                                                                                          
  void get_dWBar_dt(bool *, GeoState &, Vec<double> &, SVec<double,dim> &, SVec<double,dim> &,
                    SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &);
                                                                                                                          
  double getTimeNorm()  {  return dt.norm(); }

// Included
  template<class Scalar, int neq>
  void addToH2(bool *, VarFcn *, Vec<double> &, SVec<double,dim> &, GenMat<Scalar,neq> &);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <TimeState.C>
#endif

#endif
