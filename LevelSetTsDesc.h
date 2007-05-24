#ifndef _LEVELSET_TS_DESC_H_
#define _LEVELSET_TS_DESC_H_

#include <TsDesc.h>

#include <IoData.h>
#include <Domain.h>
#include <LevelSet.h>

struct DistInfo;

class GeoSource;
template<class Scalar, int dim> class DistSVec;

                                                                                                                

//------------------------------------------------------------------------

template<int dim>
class LevelSetTsDesc : public TsDesc<dim> {

 protected:
  LevelSet *LS;
  DistVec<double> Phi; //conservative variables
  DistVec<double> PhiV;//primitive variables
  DistSVec<double,dim> Vg;
  DistSVec<double,dim> *Vgf;

	// frequency for reinitialization of level set
	int frequencyLS;
  
 public:
  LevelSetTsDesc(IoData &, GeoSource &, Domain *);
  ~LevelSetTsDesc();

  
  LevelSet *getLevelSet() {return LS;};
  
  //-- overrides the functions implemented in TsDesc.
  void setupTimeStepping(DistSVec<double,dim> *, IoData &);
  double computeTimeStep(int, double *, DistSVec<double,dim> &);
  void updateStateVectors(DistSVec<double,dim> &, int);
  int checkSolution(DistSVec<double,dim> &);
  void setupOutputToDisk(IoData &, bool *, int, double, 
                        DistSVec<double,dim> &);
  void outputToDisk(IoData &, bool*, int, int, int, double, double, 
		    	DistSVec<double,dim> &);


  virtual int solveNonLinearSystem(DistSVec<double,dim> &)=0;

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <LevelSetTsDesc.C>
#endif

#endif
