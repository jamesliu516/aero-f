#ifndef _CURVATURE_DETECTION_H_
#define _CURVATURE_DETECTION_H_

class Domain;
class SubDomain;
template<class T> class CommPattern;
template<class Scalar> class DistVec;
template<class Scalar, int dim> class DistSVec;

//------------------------------------------------------------------------------

class CurvatureDetection {

  int numLocSub;
  SubDomain** subDomain;
//  DistVec<double>* tag;
  DistSVec<double,6>* normals;
  DistSVec<double,5>* originAndTag;
//  CommPattern<double>* vec1;
  CommPattern<double>* vec5;
  CommPattern<double>* vec6;

public:
  
  CurvatureDetection(Domain*);
  ~CurvatureDetection();

  void compute(double, int, double,  DistSVec<double,3>&, DistVec<bool>&);

};

//------------------------------------------------------------------------------

#endif
