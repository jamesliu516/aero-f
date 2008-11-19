#include <EulerStructGhostFluid.h>
#include <stdio.h>
//---------------------------------------------------------------------

EulerStructGhostFluid::EulerStructGhostFluid(SubDomain* subD, Vec<double> & phi_input, SVec<double,3> &gradphi_input, Vec<bool> &realNodeTag0_input, Vec<bool> &realNodeTag_input): philevel(phi_input), gradPhilevel(gradphi_input), realNodeTag0(realNodeTag0_input), realNodeTag(realNodeTag_input)
{
  subDomain = subD;
  triaSurf = new TriangulatedSurface;

  //phi = new double[subDomain->numNodes()];
  n_phi = new Vec3D[subDomain->numNodes()];
//  gNodeList = new int[subDomain->numNodes()];

  for (int iNode=0; iNode<(subDomain->numNodes()); iNode++)
  {
    n_phi[iNode] = 0.0;
  }

/*  tagTets = new int[subDomain->numElems()];
  for (int iElem=0; iElem<(subDomain->numElems()); iElem++)
    tagTets[iElem] = 1;
*/
//  numGNode = 0;
//  imageList = 0;
//  tetList = 0;
//  localCoord = 0;

//  tagNodes = 0;

  PhysBAM_Interface = 0;

  nodeTag = 0;
  compressible_fluid_particle = 0;
  tetrahedralized_volume = 0;
  numChosenNodes = 0;
  tempNodeList = 0;
  numChosenElems = 0;
  tempElemList = 0;
  chainMail = 0;
  solids_particle = 0;
  elem_list = 0;
  tetrahedron_mesh = 0;

//  CFPtoGlobNodeMap = 0;  
}

//---------------------------------------------------------------------------------

EulerStructGhostFluid::~EulerStructGhostFluid()
{
  if (subDomain) delete[] subDomain;
//  if (phi) delete[] phi;
  if (n_phi) delete[] n_phi;
//  if (gNodeList) delete[] gNodeList;
//  if (imageList) delete[] imageList;
//  if (tetList) delete[] tetList;
//  if (localCoord) delete[] localCoord;
//  if (tagTets) delete[] tagTets;
  if (triaSurf) delete triaSurf;
//  if (tagNodes) delete[] tagNodes;
  if (PhysBAM_Interface) delete PhysBAM_Interface;

  if (nodeTag) delete[] nodeTag;
  if (tempNodeList) delete[] tempNodeList;
  if (tempElemList) delete[] tempElemList;
  if (chainMail) delete[] chainMail;
  if (compressible_fluid_particle) delete compressible_fluid_particle;
  if (tetrahedralized_volume) delete tetrahedralized_volume;
  if (solids_particle) delete solids_particle;
  if (elem_list) delete elem_list;
  if (tetrahedron_mesh) delete tetrahedron_mesh;

}

//-----------------------------------------------------------------------------------------

void EulerStructGhostFluid::getTriangulatedSurfaceFromFace(SVec<double,3> &X)
{
  subDomain->getTriangulatedSurfaceFromFace(X,triaSurf);
}

//----------------------------------------------------------------------------------------

void EulerStructGhostFluid::getTriangulatedSurfaceFromFace()
{
  subDomain->getTriangulatedSurfaceFromFace(triaSurf);
}

//-----------------------------------------------------------------------------------------

void EulerStructGhostFluid::getPhiFromModule(SVec<double,3> &X, bool updateRealNodeTag = false)
{
  if (updateRealNodeTag) {
    for (int i=0; i<X.size(); i++) {
      realNodeTag0[i] = realNodeTag[i];
      realNodeTag[i] = true;
    }
  }

  for (int i=0; i<X.size(); i++){
    PhysBAM::VECTOR<double,3> position(X[i][0], X[i][1], X[i][2]);
    philevel[i] = PhysBAM_Interface->Phi(position);
    if (updateRealNodeTag && philevel[i]<0.0)  realNodeTag[i] = false;
  }
}

//----------------------------------------------------------------------------------------
void EulerStructGhostFluid::getPhiFromModule(SVec<double,3> &X, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, bool updateRealNodeTag = false)
{
  if (updateRealNodeTag) {
    for (int i=0; i<X.size(); i++) {
      realNodeTag0[i] = realNodeTag[i];
      realNodeTag[i] = true;
    }
  }

  for (int i=0; i<X.size(); i++){
//    double temp = philevel[i];
    philevel[i] = 10.0;
    if (X[i][0]<xmin || X[i][0]>xmax || X[i][1]<ymin || X[i][1]>ymax || X[i][2]<zmin || X[i][2]>zmax)
      continue;
    PhysBAM::VECTOR<double,3> position(X[i][0], X[i][1], X[i][2]);
//    fprintf(stderr,"node %d (%f, %f, %f).\n", i, X[i][0], X[i][1], X[i][2]);
    philevel[i] = PhysBAM_Interface->Phi(position);
//    fprintf(stderr,"old ls = %f, new ls = %f.\n", temp, philevel[i]);
    if (updateRealNodeTag && philevel[i]<0.0)  realNodeTag[i] = false;
  }
}


//-----------------------------------------------------------------------------------------

void EulerStructGhostFluid::getGradPhiFromModule(SVec<double,3> &X)
{
  for (int i=0; i<X.size(); i++){
    PhysBAM::VECTOR<double,3> position(X[i][0], X[i][1], X[i][2]);
    PhysBAM::VECTOR<double,3> gradPhi = PhysBAM_Interface->Gradient(position);
    n_phi[i][0] = gradPhi[1];  n_phi[i][1] = gradPhi[2];  n_phi[i][2] = gradPhi[3];
    gradPhilevel[i][0] = gradPhi[1]; gradPhilevel[i][1] = gradPhi[2]; gradPhilevel[i][2] = gradPhi[3];
  }
}

//-----------------------------------------------------------------------------------------

void EulerStructGhostFluid::getGradPhiFromModule(SVec<double,3> &X, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
  for (int i=0; i<X.size(); i++){
    gradPhilevel[i][0] = gradPhilevel[i][1] = gradPhilevel[i][2] = 0.0;

    if (X[i][0]<xmin || X[i][0]>xmax || X[i][1]<ymin || X[i][1]>ymax || X[i][2]<zmin || X[i][2]>zmax)
      continue;

    PhysBAM::VECTOR<double,3> position(X[i][0], X[i][1], X[i][2]);
    PhysBAM::VECTOR<double,3> gradPhi = PhysBAM_Interface->Gradient(position);
    n_phi[i][0] = gradPhi[1];  n_phi[i][1] = gradPhi[2];  n_phi[i][2] = gradPhi[3];
    gradPhilevel[i][0] = gradPhi[1]; gradPhilevel[i][1] = gradPhi[2]; gradPhilevel[i][2] = gradPhi[3];
  }
}

//----------------------------------------------------------------------------------------

Vec3D EulerStructGhostFluid::getGradPhiAtPosition(Vec3D pos)
{
  PhysBAM::VECTOR<double,3> position(pos[0], pos[1], pos[2]);
  PhysBAM::VECTOR<double,3> gradPhi = PhysBAM_Interface->Gradient(position);
  Vec3D grad(gradPhi[1], gradPhi[2], gradPhi[3]);
  return grad;
}

//-----------------------------------------------------------------------------------------

/*void EulerStructGhostFluid::getDistManually(SVec<double,3> &X)
{

  for (int i=0; i<X.size(); i++){
    phi[i] = X[i][0] - 0.6;
    n_phi[i][0] = 1.0;  n_phi[i][1] = 0.0;  n_phi[i][2] = 0.0;
  }
  
}*/
//-----------------------------------------------------------------------------------------
/*
void EulerStructGhostFluid::getListOfGhostNodes()
{
  subDomain->getGhostNodes(philevel.v, gNodeList, numGNode);

}
*/
//-----------------------------------------------------------------------------------------
/*
void EulerStructGhostFluid::getImageList(SVec<double,3> &X)
{
  if (imageList) delete[] imageList;
  imageList = new Vec3D[numGNode];
  for (int i=0; i<numGNode; i++){
    imageList[i][0] = X[gNodeList[i]][0] - 2.0*philevel[gNodeList[i]]*n_phi[gNodeList[i]][0];
    imageList[i][1] = X[gNodeList[i]][1] - 2.0*philevel[gNodeList[i]]*n_phi[gNodeList[i]][1];
    imageList[i][2] = X[gNodeList[i]][2] - 2.0*philevel[gNodeList[i]]*n_phi[gNodeList[i]][2];
  }
  fprintf(stderr, "numGNodes = %d\n", numGNode);
  for (int i=0; i<numGNode; i++)
    fprintf(stderr,"imageList[%d] = [%.4f, %.4f, %.4f]\n",i,imageList[i][0],imageList[i][1],imageList[i][2]);
}
*/
//-----------------------------------------------------------------------------------------
/*
void EulerStructGhostFluid::getTetList(SVec<double,3> &X)
{
  if (tetList) delete[] tetList;
  tetList = new int[numGNode];
  for (int iNode=0; iNode<numGNode; iNode++){

    for (int iElem=0; iElem<(subDomain->numElems()); iElem++)
      tagTets[iElem] = 1;

    int depth = 1;
    bool found = false;
 
    while (depth) {
      int listSize = 0;
      int* list = subDomain->getNeiElemOfNode(gNodeList[iNode], depth, listSize);
      fprintf(stderr,"listSize = %d\n",listSize);
      for (int i=0; i<listSize; i++){
        if (tagTets[list[i]] == 0) continue;
        found = subDomain->isINodeinITet(imageList[iNode], list[i], X);
        tagTets[list[i]] = 0;
        if (found == true) {tetList[iNode] = list[i];  break; }
      }
      if (found == true) break;      
     
      bool complete = true;
      for (int l=0; l<(subDomain->numElems()); l++)
        if (tagTets[l] == 1) {complete = false; break; }
      if (complete == true){
        fprintf(stderr, "gone over all the elements but failed..........");
        break;
      }
      
      depth++;
    }
  }
  
//  fprintf(stderr, "numGNode = %d\n", numGNode);
//  for (int i=0; i<numGNode; i++)
//    fprintf(stderr, "tetList[%d] = %d\n", i, tetList[i]);
}
*/
//----------------------------------------------------------------------------------------
/*
void EulerStructGhostFluid::getLocalCoord(SVec<double,3> &X)
{
  if (localCoord) delete[] localCoord;
  localCoord = new Vec3D[numGNode];
  for (int i=0; i<numGNode; i++)
    subDomain->localCoord(imageList[i], tetList[i], X, localCoord[i]);

  //[for debug]
//  for (int i=0; i<numGNode; i++)
//    fprintf(stderr, "gNodeList[%d]: local coordinates of its image in element [%d]: [%.4f, %.4f, %.4f ].\n", i, tetList[i], localCoord[i][0], localCoord[i][1], localCoord[i][2]);
}
*/
//---------------------------------------------------------------------------------------
/*
void EulerStructGhostFluid::getPhiValuesToWrite(SVec<double,1> &phiWrite)
{
  for (int i=0; i<subDomain->numNodes(); i++)  
    phiWrite[i][0] = philevel[i];
}
*/
//----------------------------------------------------------------------------------------
/*
void EulerStructGhostFluid::solve3by3LinearSystem(double a11, double a12, double a13,
                                                  double a21, double a22, double a23,
                                                  double a31, double a32, double a33,
                                                  double &x1, double &x2, double &x3,
                                                  double b1,  double b2,  double b3)
{
  double c11, c12, c13, c21, c22, c23, c31, c32, c33;
  double detA;      //info. for solving a 3 by 3 linear system.

  detA = a11*(a22*a33 - a32*a23) - a21*(a12*a33 - a32*a13) + a31*(a12*a23 - a22*a13);
  if (detA == 0.0) {fprintf(stderr,"error message: det(A) = 0 in EulerStructGhostFluid::solve3by3LinearSystem(...) \n");  exit(1);}

  c11 = 1.0/detA*(a22*a33-a32*a23);  c21 = 1.0/detA*(a23*a31-a33*a21);  c31 = 1.0/detA*(a21*a32-a31*a22);
  c12 = 1.0/detA*(a13*a32-a33*a12);  c22 = 1.0/detA*(a11*a33-a31*a13);  c32 = 1.0/detA*(a12*a31-a32*a11);
  c13 = 1.0/detA*(a12*a23-a22*a13);  c23 = 1.0/detA*(a13*a21-a23*a11);  c33 = 1.0/detA*(a11*a22-a12*a21);

  x1 = c11*b1 + c12*b2 + c13*b3;
  x2 = c21*b1 + c22*b2 + c23*b3;
  x3 = c31*b1 + c32*b2 + c33*b3;

}
*/
//----------------------------------------------------------------------------------------

int (*EulerStructGhostFluid::getTriangleNodeNumFromTriaSurf() const)[3]
{
  return triaSurf->getTriangleNodeNum();
}

//-----------------------------------------------------------------------------------------

int EulerStructGhostFluid::getNumTriangleFromTriaSurf()
{
  return triaSurf->getNumTriangle();
}

//------------------------------------------------------------------------------------------
/*
double EulerStructGhostFluid::tagGhostTets() //send back phi_max.
{
  //1. initialize a tag on nodes.
  if (!tagNodes) tagNodes = new int[subDomain->numNodes()];
  for (int iNode=0; iNode<(subDomain->numNodes()); iNode++)
    tagNodes[iNode] = 0;

  //2. tag the "ghost nodes".
  getListOfGhostNodes();
  for (int i=0; i<numGNode; i++)
    tagNodes[gNodeList[i]] = 1;

  //3. tag all neighbour elements of ghost nodes. find phiMax.
  if (tagTets)
    for (int iElem=0; iElem<(subDomain->numElems()); iElem++)
      tagTets[iElem] = 0;
  if (!tagTets){
    tagTets = new int[subDomain->numElems()];
    for (int iElem=0; iElem<(subDomain->numElems()); iElem++)
    tagTets[iElem] = 0;
  }

  double phiMax = 0.0;
  for (int iElem=0; iElem<(subDomain->numElems()); iElem++) {
    int* elemNodes = subDomain->getElemNodeNum(iElem);
    for (int j=0; j<4; j++)
      if (tagNodes[elemNodes[j]] == 1) 
        tagTets[iElem] = 1;  //tag the neighbour elements of ghost nodes.

    if (tagTets[iElem] == 1) //update phiMax.
      for (int j=0; j<4; j++)
        if (philevel[elemNodes[j]] > phiMax)  phiMax = philevel[elemNodes[j]];
  }
    
  return phiMax;

}
*/
//---------------------------------------------------------------------------------------
/*
void EulerStructGhostFluid::tagRealFluidTets(double phimax)
{
  for (int iElem=0; iElem<(subDomain->numElems()); iElem++) {
    int* elemNodes = subDomain->getElemNodeNum(iElem);
    if (   (philevel[elemNodes[0]]>0 && philevel[elemNodes[1]]>0 && philevel[elemNodes[2]]>0 )
        && (philevel[elemNodes[0]]<phimax || philevel[elemNodes[1]]<phimax || philevel[elemNodes[2]]<phimax)  )
      tagTets[iElem] = 1;
  }
}
*/
//----------------------------------------------------------------------------------------
/*
int EulerStructGhostFluid::reconstructTagNodes()
{
  int num_nodes_input = 0;
  if (!tagNodes) tagNodes = new int[subDomain->numNodes()];
  for (int iNode=0; iNode<(subDomain->numNodes()); iNode++)
    tagNodes[iNode] = -1;
  
  for (int iElem=0; iElem<(subDomain->numElems()); iElem++){
    if (tagTets[iElem]==1) {
      int* elemNodes = subDomain->getElemNodeNum(iElem);
      for (int j=0; j<4; j++)
        if (tagNodes[elemNodes[j]]<0) {tagNodes[elemNodes[j]] = num_nodes_input;  num_nodes_input++; }
    }
  }
  return num_nodes_input;
}
*/
//-------------------------------------------------------------------------------------------
/*
double EulerStructGhostFluid::phiAtNode(const int iNode)
{
  return philevel[iNode];
}
*/
//--------------------------------------------------------------------------------------------

double EulerStructGhostFluid::specifyBandwidth()
{
  return subDomain->specifyBandwidth(philevel);
}

//-------------------------------------------------------------------------------------------

void EulerStructGhostFluid::constructInterface(PhysBAM::GRID_3D<double> grid3d) 
{
  if (PhysBAM_Interface) delete PhysBAM_Interface;
  PhysBAM_Interface = new PhysBAM::AERO_INTERFACE_1<double> (grid3d);
}

//--------------------------------------------------------------------------------------------

void EulerStructGhostFluid::computeLevelSet(PhysBAM::TRIANGULATED_SURFACE<double> &physbam_triangulated_surface)
{
  PhysBAM_Interface->Compute_Level_Set(physbam_triangulated_surface);
}

//--------------------------------------------------------------------------------------------

bool EulerStructGhostFluid::insideOutside(double *position, const double xmin, const double xmax,
                     const double ymin, const double ymax, const double zmin, const double zmax)
{
  if (position[0] < xmin || position[0] > xmax || position[1] < ymin || position[1] > ymax ||
      position[2] < zmin || position[2] > zmax )  return false;
  return true;
}

//---------------------------------------------------------------------------------------------
bool EulerStructGhostFluid::insideOutside(Vec3D position, const double xmin, const double xmax,
                     const double ymin, const double ymax, const double zmin, const double zmax)
{
  if (position[0] < xmin || position[0] > xmax || position[1] < ymin || position[1] > ymax ||
      position[2] < zmin || position[2] > zmax )  return false;
  return true;
}

//-------------------------------------------------------------------------------

void EulerStructGhostFluid::initializePhysBAMMPI(int numNodes,double depth, int cpuNum)
{
  int* locToGlob = subDomain->getNodeMap();
  PhysBAM::ARRAY<int> CFPtoGlobNodeMap(numChosenNodes);
  for (int iParticle=1; iParticle<=numChosenNodes; iParticle++)
    CFPtoGlobNodeMap(iParticle)=locToGlob[tempNodeList[iParticle-1]]+1;

  PhysBAM::VECTOR<double,3> dimCartMeshMin, dimCartMeshMax; 
  PhysBAM_Interface->Initialize_MPI(*tetrahedralized_volume,CFPtoGlobNodeMap,numNodes,depth, dimCartMeshMin, dimCartMeshMax); 

}

//-------------------------------------------------------------------------------

void EulerStructGhostFluid::computeGhostNodes()
{
  if ((!PhysBAM_Interface)||(!solids_particle)||(!tetrahedralized_volume)||(!elem_list)||(!tetrahedron_mesh)) {
    fprintf(stderr,"ingredients lost!\n"); return;
  }

  PhysBAM_Interface->Compute_Ghost_Cells(*tetrahedralized_volume, *compressible_fluid_particle);

}

//-------------------------------------------------------------------------------

void EulerStructGhostFluid::updateRealNodeTag(SVec<double,3>& X, double xmin, double xmax, double ymin,
                                              double ymax, double zmin, double zmax)
{
  for (int i=0; i<subDomain->numNodes(); i++) {  
    if (X[i][0]<xmin || X[i][0]>xmax || X[i][1]<ymin || X[i][1]>ymax || X[i][2]<zmin || X[i][2]>zmax)
      continue;
    realNodeTag0[i] = realNodeTag[i];
    if (philevel[i]>=0.0)  realNodeTag[i] = true;  else realNodeTag[i] = false;
  }
}

//-------------------------------------------------------------------------------

void EulerStructGhostFluid::recomputeLevelSet(Vec3D dsp)
{
  PhysBAM::VECTOR<double,3> physbam_dsp(dsp[0],dsp[1],dsp[2]);
//  fprintf(stderr,"physbam displacement: %f, %f, %f.\n", physbam_dsp[1], physbam_dsp[2], physbam_dsp[3]); 
  PhysBAM_Interface->Frame().t += physbam_dsp; 
}





















