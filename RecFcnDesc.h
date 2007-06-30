#ifndef _REC_FCN_DESC_H_
#define _REC_FCN_DESC_H_

#include <RecFcn.h>

//------------------------------------------------------------------------------

template<int dim>
class RecFcnConstant : public RecFcn {

public:

  RecFcnConstant() : RecFcn(0.0, 0.0) {}
  ~RecFcnConstant() {}

  void precompute(double *, double *, double *, double *, 
		  double *, double *, double *, double *);
  void compute(double *, double *, double *, double *, double *, double *);


};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnLinear : public RecFcn {

public:

  RecFcnLinear(double b, double e) : RecFcn(b, e) {}
  ~RecFcnLinear() {}

  void precompute(double *, double *, double *, double *, 
		  double *, double *, double *, double *); 
  void compute(double *, double *, double *, double *, double *, double *);

};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnVanAlbada : public RecFcn {

public:

  RecFcnVanAlbada(double b, double e) : RecFcn(b, e) {}
  ~RecFcnVanAlbada() {}

  void precompute(double *, double *, double *, double *, 
		  double *, double *, double *, double *); 
  void compute(double *, double *, double *, double *, double *, double *);

};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnLtdMultiDim : public RecFcnLinear<dim> {

public:

  RecFcnLtdMultiDim(double b, double e) : RecFcnLinear<dim>(b, e) {}
  ~RecFcnLtdMultiDim() {}

  virtual void computeLimiter(double *, double *, double *, double *, double,
			      double *, double *, double *, double *, double,
			      double *, double *) = 0;

};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnBarth : public RecFcnLtdMultiDim<dim> {

public:

  RecFcnBarth(double b, double e) : RecFcnLtdMultiDim<dim>(b, e) {}
  ~RecFcnBarth() {}

  void computeLimiter(double *, double *, double *, double *, double,
		      double *, double *, double *, double *, double,
		      double *, double *);
};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnVenkat : public RecFcnLtdMultiDim<dim> {

public:

  RecFcnVenkat(double b, double e) : RecFcnLtdMultiDim<dim>(b, e) {}
  ~RecFcnVenkat() {}

  void computeLimiter(double *, double *, double *, double *, double,
		      double *, double *, double *, double *, double,
		      double *, double *);

};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnLinearConstant : public RecFcn {

public:

  RecFcnLinearConstant(double b, double e) : RecFcn(b, e) {}
  ~RecFcnLinearConstant() {}

  void compute(double *, double *, double *, double *, double *, double *);

};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnVanAlbadaConstant : public RecFcn {

public:

  RecFcnVanAlbadaConstant(double b, double e) : RecFcn(b, e) {}
  ~RecFcnVanAlbadaConstant() {}

  void compute(double *, double *, double *, double *, double *, double *);

};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnLinearVanAlbada : public RecFcn {

public:

  RecFcnLinearVanAlbada(double b, double e) : RecFcn(b, e) {}
  ~RecFcnLinearVanAlbada() {}

  void compute(double *, double *, double *, double *, double *, double *);

};

//------------------------------------------------------------------------------

class RecFcnLtdSensor {

  double threshold;

public:

  RecFcnLtdSensor(double eps) { threshold = eps; }
  ~RecFcnLtdSensor() {}
  
  virtual double getThreshold() { return threshold; }

};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnLtdLinear : public RecFcnLinear<dim>, public RecFcnLtdSensor {

public:

  RecFcnLtdLinear(double b, double e) : RecFcnLinear<dim>(b, e), RecFcnLtdSensor(e) {}
  ~RecFcnLtdLinear() {}

};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnLtdLinearConstant : public RecFcnLinearConstant<dim>, public RecFcnLtdSensor {

public:

  RecFcnLtdLinearConstant(double b, double e) : RecFcnLinearConstant<dim>(b, e), RecFcnLtdSensor(e) {}
  ~RecFcnLtdLinearConstant() {}

};

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnConstant<dim>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				     double* aij, double* aji, double* bij, double* bji)
{

  for (int k=0; k<dim; ++k)
    preconstant(aij[k], aji[k], bij[k], bji[k]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnConstant<5>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				   double* aij, double* aji, double* bij, double* bji)
{

  preconstant(aij[0], aji[0], bij[0], bji[0]);
  preconstant(aij[1], aji[1], bij[1], bji[1]);
  preconstant(aij[2], aji[2], bij[2], bji[2]);
  preconstant(aij[3], aji[3], bij[3], bji[3]);
  preconstant(aij[4], aji[4], bij[4], bji[4]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnConstant<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				  double* Vij, double* Vji)
{

  for (int k=0; k<dim; ++k)
    constant(Vi[k], Vj[k], Vij[k], Vji[k]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnConstant<5>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				double* Vij, double* Vji)
{

  constant(Vi[0], Vj[0], Vij[0], Vji[0]);
  constant(Vi[1], Vj[1], Vij[1], Vji[1]);
  constant(Vi[2], Vj[2], Vij[2], Vji[2]);
  constant(Vi[3], Vj[3], Vij[3], Vji[3]);
  constant(Vi[4], Vj[4], Vij[4], Vji[4]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnConstant<6>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				double* Vij, double* Vji)
{

  constant(Vi[0], Vj[0], Vij[0], Vji[0]);
  constant(Vi[1], Vj[1], Vij[1], Vji[1]);
  constant(Vi[2], Vj[2], Vij[2], Vji[2]);
  constant(Vi[3], Vj[3], Vij[3], Vji[3]);
  constant(Vi[4], Vj[4], Vij[4], Vji[4]);
  constant(Vi[5], Vj[5], Vij[5], Vji[5]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnConstant<7>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				double* Vij, double* Vji)
{

  constant(Vi[0], Vj[0], Vij[0], Vji[0]);
  constant(Vi[1], Vj[1], Vij[1], Vji[1]);
  constant(Vi[2], Vj[2], Vij[2], Vji[2]);
  constant(Vi[3], Vj[3], Vij[3], Vji[3]);
  constant(Vi[4], Vj[4], Vij[4], Vji[4]);
  constant(Vi[5], Vj[5], Vij[5], Vji[5]);
  constant(Vi[6], Vj[6], Vij[6], Vji[6]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnLinear<dim>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				   double* aij, double* aji, double* bij, double* bji)
{

  for (int k=0; k<dim; ++k)
    prelinear(aij[k], aji[k], bij[k], bji[k]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnLinear<5>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				 double* aij, double* aji, double* bij, double* bji)
{

  prelinear(aij[0], aji[0], bij[0], bji[0]);
  prelinear(aij[1], aji[1], bij[1], bji[1]);
  prelinear(aij[2], aji[2], bij[2], bji[2]);
  prelinear(aij[3], aji[3], bij[3], bji[3]);
  prelinear(aij[4], aji[4], bij[4], bji[4]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnLinear<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				double* Vij, double* Vji)
{

  for (int k=0; k<dim; ++k)
    linear(Vi[k], ddVij[k], Vj[k], ddVji[k], Vij[k], Vji[k]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnLinear<5>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
			      double* Vij, double* Vji)
{

  linear(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  linear(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  linear(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  linear(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  linear(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnVanAlbada<dim>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				      double* aij, double* aji, double* bij, double* bji)
{

  for (int k=0; k<dim; ++k)
    prevanalbada(Vi[k], ddVij[k], Vj[k], ddVji[k], aij[k], aji[k], bij[k], bji[k]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnVanAlbada<5>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				    double* aij, double* aji, double* bij, double* bji)
{

  prevanalbada(Vi[0], ddVij[0], Vj[0], ddVji[0], aij[0], aji[0], bij[0], bji[0]);
  prevanalbada(Vi[1], ddVij[1], Vj[1], ddVji[1], aij[1], aji[1], bij[1], bji[1]);
  prevanalbada(Vi[2], ddVij[2], Vj[2], ddVji[2], aij[2], aji[2], bij[2], bji[2]);
  prevanalbada(Vi[3], ddVij[3], Vj[3], ddVji[3], aij[3], aji[3], bij[3], bji[3]);
  prevanalbada(Vi[4], ddVij[4], Vj[4], ddVji[4], aij[4], aji[4], bij[4], bji[4]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnVanAlbada<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				   double* Vij, double* Vji)
{

  for (int k=0; k<dim; ++k)
    vanalbada(Vi[k], ddVij[k], Vj[k], ddVji[k], Vij[k], Vji[k]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnVanAlbada<5>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				 double* Vij, double* Vji)
{

  vanalbada(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  vanalbada(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  vanalbada(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  vanalbada(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  vanalbada(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnVanAlbada<6>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				 double* Vij, double* Vji)
{

  vanalbada(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  vanalbada(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  vanalbada(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  vanalbada(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  vanalbada(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
  vanalbada(Vi[5], ddVij[5], Vj[5], ddVji[5], Vij[5], Vji[5]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnVanAlbada<7>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				 double* Vij, double* Vji)
{

  vanalbada(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  vanalbada(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  vanalbada(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  vanalbada(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  vanalbada(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
  vanalbada(Vi[5], ddVij[5], Vj[5], ddVji[5], Vij[5], Vji[5]);
  vanalbada(Vi[6], ddVij[6], Vj[6], ddVji[6], Vij[6], Vji[6]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnBarth<dim>::computeLimiter(double *Vimax, double *Vimin, double *Vi, 
				      double *Vij, double ctrlVoli,
				      double *Vjmax, double *Vjmin, double *Vj, 
				      double *Vji, double ctrlVolj,
				      double *phii, double *phij)
{

  for (int k=0; k<dim; ++k)
    this->barth(Vimax[k], Vimin[k], Vi[k], Vij[k], ctrlVoli,
	  Vjmax[k], Vjmin[k], Vj[k], Vji[k], ctrlVolj, phii[k], phij[k]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnBarth<5>::computeLimiter(double *Vimax, double *Vimin, double *Vi, 
				    double *Vij, double ctrlVoli,
				    double *Vjmax, double *Vjmin, double *Vj, 
				    double *Vji, double ctrlVolj,
				    double *phii, double *phij)
{

  this->barth(Vimax[0], Vimin[0], Vi[0], Vij[0], ctrlVoli,
	Vjmax[0], Vjmin[0], Vj[0], Vji[0], ctrlVolj, phii[0], phij[0]);
  this->barth(Vimax[1], Vimin[1], Vi[1], Vij[1], ctrlVoli,
	Vjmax[1], Vjmin[1], Vj[1], Vji[1], ctrlVolj, phii[1], phij[1]);
  this->barth(Vimax[2], Vimin[2], Vi[2], Vij[2], ctrlVoli,
	Vjmax[2], Vjmin[2], Vj[2], Vji[2], ctrlVolj, phii[2], phij[2]);
  this->barth(Vimax[3], Vimin[3], Vi[3], Vij[3], ctrlVoli,
	Vjmax[3], Vjmin[3], Vj[3], Vji[3], ctrlVolj, phii[3], phij[3]);
  this->barth(Vimax[4], Vimin[4], Vi[4], Vij[4], ctrlVoli,
	Vjmax[4], Vjmin[4], Vj[4], Vji[4], ctrlVolj, phii[4], phij[4]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnVenkat<dim>::computeLimiter(double *Vimax, double *Vimin, double *Vi, 
				       double *Vij, double ctrlVoli,
				       double *Vjmax, double *Vjmin, double *Vj, 
				       double *Vji, double ctrlVolj,
				       double *phii, double *phij)
{

  for (int k=0; k<dim; ++k)
    this->venkat(Vimax[k], Vimin[k], Vi[k], Vij[k], ctrlVoli,
	   Vjmax[k], Vjmin[k], Vj[k], Vji[k], ctrlVolj, phii[k], phij[k]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnVenkat<5>::computeLimiter(double *Vimax, double *Vimin, double *Vi, 
				     double *Vij, double ctrlVoli,
				     double *Vjmax, double *Vjmin, double *Vj, 
				     double *Vji, double ctrlVolj,
				     double *phii, double *phij)
{

  this->venkat(Vimax[0], Vimin[0], Vi[0], Vij[0], ctrlVoli,
	 Vjmax[0], Vjmin[0], Vj[0], Vji[0], ctrlVolj, phii[0], phij[0]);
  this->venkat(Vimax[1], Vimin[1], Vi[1], Vij[1], ctrlVoli,
	 Vjmax[1], Vjmin[1], Vj[1], Vji[1], ctrlVolj, phii[1], phij[1]);
  this->venkat(Vimax[2], Vimin[2], Vi[2], Vij[2], ctrlVoli,
	 Vjmax[2], Vjmin[2], Vj[2], Vji[2], ctrlVolj, phii[2], phij[2]);
  this->venkat(Vimax[3], Vimin[3], Vi[3], Vij[3], ctrlVoli,
	 Vjmax[3], Vjmin[3], Vj[3], Vji[3], ctrlVolj, phii[3], phij[3]);
  this->venkat(Vimax[4], Vimin[4], Vi[4], Vij[4], ctrlVoli,
	 Vjmax[4], Vjmin[4], Vj[4], Vji[4], ctrlVolj, phii[4], phij[4]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnLinearConstant<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
					double* Vij, double* Vji)
{

  fprintf(stderr, "*** Error: RecFcnLinearConstant<%d>::compute is not overloaded\n", dim);
  exit(1);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnLinearConstant<6>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				      double* Vij, double* Vji)
{

  linear(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  linear(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  linear(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  linear(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  linear(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
  constant(Vi[5], Vj[5], Vij[5], Vji[5]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnLinearConstant<7>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				      double* Vij, double* Vji)
{

  linear(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  linear(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  linear(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  linear(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  linear(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
  constant(Vi[5], Vj[5], Vij[5], Vji[5]);
  constant(Vi[6], Vj[6], Vij[6], Vji[6]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnVanAlbadaConstant<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
					   double* Vij, double* Vji)
{

  fprintf(stderr, "*** Error: RecFcnVanAlbadaConstant<%d>::compute is not overloaded\n", dim);
  exit(1);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnVanAlbadaConstant<6>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
					 double* Vij, double* Vji)
{

  vanalbada(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  vanalbada(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  vanalbada(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  vanalbada(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  vanalbada(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
  constant(Vi[5], Vj[5], Vij[5], Vji[5]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnVanAlbadaConstant<7>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
					 double* Vij, double* Vji)
{

  vanalbada(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  vanalbada(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  vanalbada(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  vanalbada(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  vanalbada(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
  constant(Vi[5], Vj[5], Vij[5], Vji[5]);
  constant(Vi[6], Vj[6], Vij[6], Vji[6]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnLinearVanAlbada<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
					 double* Vij, double* Vji)
{

  fprintf(stderr, "*** Error: RecFcnLinearVanAlbada<%d>::compute is not overloaded\n", dim);
  exit(1);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnLinearVanAlbada<6>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				       double* Vij, double* Vji)
{

  linear(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  linear(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  linear(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  linear(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  linear(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
  vanalbada(Vi[5], ddVij[5], Vj[5], ddVji[5], Vij[5], Vji[5]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnLinearVanAlbada<7>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				       double* Vij, double* Vji)
{

  linear(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  linear(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  linear(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  linear(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  linear(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
  vanalbada(Vi[5], ddVij[5], Vj[5], ddVji[5], Vij[5], Vji[5]);
  vanalbada(Vi[6], ddVij[6], Vj[6], ddVji[6], Vij[6], Vji[6]);

}

//------------------------------------------------------------------------------

#endif
