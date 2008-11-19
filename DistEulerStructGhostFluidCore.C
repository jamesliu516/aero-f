#include <DistEulerStructGhostFluid.h>
#include <stdio.h>
//------------------------------------------------------------------------------------

DistEulerStructGhostFluid::DistEulerStructGhostFluid(Domain* domain, IoData& iod):numGlobNodes(domain->numNodes())
{
  numLocSub = domain->getNumLocSub();
  subDomain = domain->getSubDomain();
//  numGlobNodes = domain->numNodes();

  PhiWrite = new DistSVec<double,1>(domain->getNodeDistInfo());

  philevel = new DistVec<double> (domain->getNodeDistInfo());
  gradPhilevel = new DistSVec<double,3>(domain->getNodeDistInfo());

  realNodeTag0 = new DistVec<bool> (domain->getNodeDistInfo());
  realNodeTag = new DistVec<bool> (domain->getNodeDistInfo());


  subESGF = new EulerStructGhostFluid *[numLocSub];

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) 
    subESGF[iSub] = new EulerStructGhostFluid(subDomain[iSub], (*philevel)(iSub), (*gradPhilevel)(iSub), 
                                              (*realNodeTag0)(iSub), (*realNodeTag)(iSub));

  length_triangle_list = 0;
  length_solids_particle_list = 0;
  triangle_list = 0;
  solids_particle_list = 0;
  tag = 0;
//  PhysBAM_Interface = 0;
//  nodeTag = 0;
//  compressible_fluid_particle = 0;
//  tetrahedralized_volume = 0;
//  numChosenNodes = 0;
//  tempNodeList = 0;
//  solids_particle = 0;
//  elem_list = 0;
//  tetrahedron_mesh = 0;

  bandwidth = 0.0;

  if (strlen(iod.input.solidsurface)>0) {
    solidsurface = new char[strlen(iod.input.prefix) + strlen(iod.input.solidsurface) + 1];
    strcpy(solidsurface, iod.input.prefix);
    strcat(solidsurface, iod.input.solidsurface);
  } else solidsurface = 0;

  if (iod.input.ghostsolid.boundingBox.provide == 1 && 
      ( (iod.input.ghostsolid.boundingBox.xmax-iod.input.ghostsolid.boundingBox.xmin> 1e-12) &&
        (iod.input.ghostsolid.boundingBox.ymax-iod.input.ghostsolid.boundingBox.ymin> 1e-12) &&
        (iod.input.ghostsolid.boundingBox.zmax-iod.input.ghostsolid.boundingBox.zmin> 1e-12) ) ) {
    givenBB = true;
    Xmin = iod.input.ghostsolid.boundingBox.xmin;
    Xmax = iod.input.ghostsolid.boundingBox.xmax;
    Ymin = iod.input.ghostsolid.boundingBox.ymin;
    Ymax = iod.input.ghostsolid.boundingBox.ymax;
    Zmin = iod.input.ghostsolid.boundingBox.zmin;
    Zmax = iod.input.ghostsolid.boundingBox.zmax;
  } else {givenBB = false; Xmin = Xmax = Ymin = Ymax = Zmin = Zmax = 0.0;}

  timer = domain->getTimer();  
  com = domain->getCommunicator();

  // for Forced Motion with a constant translating velocity
  vel = 0.0;  //vel[0] = 10.0/iod.ref.rv.velocity; //TODO: in future this should be loaded from the input file.
//  vel[0] = 10.0/iod.ref.rv.velocity;
  dsp = 0.0;


  if (iod.rmesh.timestep>0) {fprintf(stderr,"rmesh exists.\n");}
  else fprintf(stderr,"rmesh doesn't exist.\n");

}

//------------------------------------------------------------------------------------

DistEulerStructGhostFluid::~DistEulerStructGhostFluid()
{
  if (subESGF) {
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; iSub++)
      if (subESGF[iSub])  delete subESGF[iSub];
    delete [] subESGF;
  }

  if (triangle_list) {
    for (int i=0; i<length_triangle_list; i++)
      if (triangle_list[i]) delete[] triangle_list[i];
    delete[] triangle_list;
  }

  if (solids_particle_list) delete[] solids_particle_list;
  

  if (tag) delete[] tag;

 // if (PhysBAM_Interface) delete PhysBAM_Interface;

  if (philevel) delete philevel;
  if (PhiWrite) delete PhiWrite;

//  if (nodeTag) delete[] nodeTag;
//  if (tempNodeList) delete[] tempNodeList;
//  if (compressible_fluid_particle) delete compressible_fluid_particle;
//  if (tetrahedralized_volume) delete tetrahedralized_volume;  
//  if (solids_particle) delete solids_particle;
//  if (elem_list) delete elem_list;
//  if (tetrahedron_mesh) delete tetrahedron_mesh;
  if (solidsurface) delete[] solidsurface;


}

//------------------------------------------------------------------------------------

void DistEulerStructGhostFluid::getBoundingBox(double &xmin, double &xmax, double &ymin, double &ymax, double &zmin, double &zmax)
{ xmin = Xmin;  xmax = Xmax;  ymin = Ymin;  ymax = Ymax;  zmin = Zmin;  zmax = Zmax; }

//-----------------------------------------------------------------------------------

void DistEulerStructGhostFluid::getTriangulatedSurfaceFromFace(DistSVec<double, 3> &X)
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->getTriangulatedSurfaceFromFace(X(iSub));
}

//------------------------------------------------------------------------------------

void DistEulerStructGhostFluid::getTriangulatedSurfaceFromFace()
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->getTriangulatedSurfaceFromFace();
}

//------------------------------------------------------------------------------------

void DistEulerStructGhostFluid::getPhiFromModule(DistSVec<double, 3> &X, bool updateRealNodeTag=false)
{
  fprintf(stderr, "DistEulerStructGhostFluid::getPhiFromModule called. updateRealNodeTag = %d.\n", updateRealNodeTag);
#pragma omp parallel for 
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->getPhiFromModule(X(iSub), updateRealNodeTag);
}

//------------------------------------------------------------------------------------

void DistEulerStructGhostFluid::getPhiFromModule(DistSVec<double, 3> &X, double xmin, double xmax,
                                                 double ymin, double ymax, double zmin, double zmax,
                                                 bool updateRealNodeTag = false)
{
  fprintf(stderr, "DistEulerStructGhostFluid::getPhiFromModule called. updateRealNodeTag = %d.\n", updateRealNodeTag);
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->getPhiFromModule(X(iSub), xmin, xmax, ymin, ymax, zmin, zmax, updateRealNodeTag);
}
//-------------------------------------------------------------------------------------

/*double DistEulerStructGhostFluid::getPhiFromModule(double x, double y, double z)
{
  PhysBAM::VECTOR<double,3> position(x, y, z);
  return PhysBAM_Interface->Phi(position);
}*/

//-------------------------------------------------------------------------------------

void DistEulerStructGhostFluid::getGradPhiFromModule(DistSVec<double, 3> &X)
{
  fprintf(stderr,"DistEulerStructGhostFluid::getGradPhiFromModule(DistSVec<double, 3> &X) called.\n");
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->getGradPhiFromModule(X(iSub));
}

//--------------------------------------------------------------------------------------

void DistEulerStructGhostFluid::getGradPhiFromModule(DistSVec<double, 3> &X, double xmin, double xmax,
                                                     double ymin, double ymax, double zmin, double zmax)
{
  fprintf(stderr,"DistEulerStructGhostFluid::getGradPhiFromModule called.\n");
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->getGradPhiFromModule(X(iSub), xmin, xmax, ymin, ymax, zmin, zmax);
}

//-----------------------------------------------------------------------------------------

/*void DistEulerStructGhostFluid::getDistManually(DistSVec<double,3> &X)
{

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->getDistManually(X(iSub));

}*/

//-------------------------------------------------------------------------------------

/*void DistEulerStructGhostFluid::getListOfGhostNodes()
{
  fprintf(stderr,"DistEulerStructGhostFluid::getListOfGhostNodes() called.\n");
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->getListOfGhostNodes();
}
*/
//--------------------------------------------------------------------------------------
/*
void DistEulerStructGhostFluid::getImageList(DistSVec<double, 3> &X)
{
  fprintf(stderr,"DistEulerStructGhostFluid::getImageList(DistSVec<double, 3> &X) called.\n");
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->getImageList(X(iSub));
}
*/
//---------------------------------------------------------------------------------------
/*
void DistEulerStructGhostFluid::getTetList(DistSVec<double, 3> &X)
{
  fprintf(stderr,"DistEulerStructGhostFluid::getTetList(DistSVec<double, 3> &X) called.\n");
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->getTetList(X(iSub));
}
*/
//----------------------------------------------------------------------------------------
/*
void DistEulerStructGhostFluid::getLocalCoord(DistSVec<double, 3> &X)
{
  fprintf(stderr,"DistEulerStructGhostFluid::getLocalCoord(DistSVec<double, 3> &X) called.\n");
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->getLocalCoord(X(iSub));
}
*/
//----------------------------------------------------------------------------------------
/*
void DistEulerStructGhostFluid::getPhiValuesToWrite()
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->getPhiValuesToWrite((*PhiWrite)(iSub));
}
*/
//----------------------------------------------------------------------------------------

void DistEulerStructGhostFluid::prepareForCommunication(DistSVec<double, 3> &X)
// construct triangle_list & solids_particle_list. Prepare for communication with PhysBAM.

{
  fprintf(stderr,"DistEulerStructGhostFluid::prepareForCommunication called.\n");
  //1. construct the tag.
  if (tag) delete[] tag;
  int length_tag = 0;
  for (int iSub=0; iSub<numLocSub; iSub++) 
    length_tag += subDomain[iSub]->numNodes();
  tag = new int[length_tag];
  for (int i=0; i<length_tag; i++)	tag[i] = -1;

  //2. delete old lists.
  if (triangle_list) {
    for (int i=0; i<length_triangle_list; i++)
      if (triangle_list[i]) delete[] triangle_list[i];
    delete[] triangle_list;
  }
  length_triangle_list = 0;

  if (solids_particle_list)  delete[] solids_particle_list;
  length_solids_particle_list = 0;

  //3. construct temporary triangle list and solids particle list.
  int length_temp_triangle_list = 0;
  int length_temp_solids_particle_list = length_tag;
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) 
    length_temp_triangle_list += subDomain[iSub]->numElems();

  int (*temp_triangle_list)[3];
  temp_triangle_list = new int[length_temp_triangle_list][3];
  double (*temp_solids_particle_list)[3];
  temp_solids_particle_list = new double[length_temp_solids_particle_list][3];

  //4. Main loop. Construct triang_list & solids_particle_list.
  for (int iSub=0; iSub<numLocSub; iSub++) {
    int numTriangle_iSub = subESGF[iSub]->getNumTriangleFromTriaSurf();
    int (*triangle_list_iSub)[3];
    triangle_list_iSub = subESGF[iSub]->getTriangleNodeNumFromTriaSurf();
    
    for (int i=0; i<numTriangle_iSub; i++) {
      for (int j=0; j<3; j++) {
        if (tag[triangle_list_iSub[i][j]] < 0 ) {
          tag[triangle_list_iSub[i][j]] = length_solids_particle_list;
          double x, y, z;
          subDomain[iSub]->getNodeCoords(triangle_list_iSub[i][j], X(iSub), x, y, z);
          temp_solids_particle_list[length_solids_particle_list][0] = x;
          temp_solids_particle_list[length_solids_particle_list][1] = y;
          temp_solids_particle_list[length_solids_particle_list][2] = z;
          length_solids_particle_list++;

          temp_triangle_list[length_triangle_list][j] = length_solids_particle_list - 1;          
        }
        else temp_triangle_list[length_triangle_list][j] = tag[triangle_list_iSub[i][j]];
      }
      length_triangle_list++;
    }
  }

  //5. Copy data from temp-lists to triangle_list and solids_particle_list;
  triangle_list = new int[length_triangle_list][3];
  for (int i=0; i<length_triangle_list; i++)
    for (int j=0; j<3; j++)
      triangle_list[i][j] = temp_triangle_list[i][j];
  
  solids_particle_list = new Vec3D[length_solids_particle_list];
  for (int i=0; i<length_solids_particle_list; i++)
    for (int j=0; j<3; j++)
      solids_particle_list[i][j] = temp_solids_particle_list[i][j];

  //6. clean up.
  if (temp_triangle_list) {
    delete[] temp_triangle_list;
  }
  
  if (temp_solids_particle_list) {
    delete[] temp_solids_particle_list;
  }

  //7, print to screen. For debug only.
  fprintf(stderr, "number of triangles: %d,  triangle_list:\n", length_triangle_list);
  for (int i=0; i<length_triangle_list; i++)
    fprintf(stderr, "Triangle # %d: [ %d, %d, %d ]. \n", i, triangle_list[i][0], triangle_list[i][1],
                     triangle_list[i][2]);
  fprintf(stderr, "\n number of solids_particles: %d,  solids_particle_list:\n", length_solids_particle_list);
  for (int i=0; i<length_solids_particle_list; i++)
    fprintf(stderr, "Solids Particle # %d: [ %.4f, %.4f, %.4f ]. \n", i, solids_particle_list[i][0],
                     solids_particle_list[i][1], solids_particle_list[i][2]);


}

//----------------------------------------------------------------------------------

void  DistEulerStructGhostFluid::prepareForCommunication()
{
  fprintf(stderr,"DistEulerStructGhostFluid::prepareForCommunication called. \n");
  //1.delete existing lists.
  if (triangle_list) {
    for (int i=0; i<length_triangle_list; i++)
      if (triangle_list[i]) delete[] triangle_list[i];
  }
  length_triangle_list = 0;

  if (solids_particle_list) delete[] solids_particle_list; 
  length_solids_particle_list = 0;

  //2.read data from "prolateSurface.top". 
  FILE *topFile;
  topFile = fopen(solidsurface, "r");
  if (topFile == NULL) {fprintf(stderr, "topFile doesn't exist at all :(\n"); exit(1); }
  fscanf(topFile,"%d %d", &length_solids_particle_list, &length_triangle_list);
  triangle_list = new int[length_triangle_list][3];
  solids_particle_list = new Vec3D[length_solids_particle_list];
  
  int thisNode;
  for (int iNode=0; iNode<length_solids_particle_list; iNode++) 
    fscanf(topFile, "%d %lf %lf %lf", &thisNode, &(solids_particle_list[iNode][0]), 
                                           &(solids_particle_list[iNode][1]), &(solids_particle_list[iNode][2]));
  if (thisNode!=length_solids_particle_list) {fprintf(stderr,"error in loading surface from file!\n"); exit(1);}

  int nothing;
  for (int iElem=0; iElem<length_triangle_list; iElem++) {
    fscanf(topFile, "%d %d %d %d %d", &thisNode, &nothing, &(triangle_list[iElem][0]), &(triangle_list[iElem][1]),
                                         &(triangle_list[iElem][2]));
    triangle_list[iElem][0]--; triangle_list[iElem][1]--; triangle_list[iElem][2]--;
  }
  if (thisNode!=length_triangle_list) {fprintf(stderr,"error in loading surface from file!\n"); exit(1);}
  fclose(topFile);

  //3. verify (1)triangulated surface is closed (2) normal's of all triangles point outward.
  fprintf(stderr,"Checking the solid surface...\n");
  if (checkTriangulatedSurface()) fprintf(stderr,"OK.\n");
  else exit(-1);
}

//------------------------------------------------------------------------------------------

bool DistEulerStructGhostFluid::checkTriangulatedSurface() {
  if ((length_triangle_list==0) || (length_solids_particle_list==0))
    {fprintf(stderr,"Solid surface not loaded.\n"); return false;}
  int numOfEdges = length_triangle_list*3/2,iEdge=0;
  if ((numOfEdges-((double)length_triangle_list)*3.0/2.0>0.1) || 
      (numOfEdges-((double)length_triangle_list)*3.0/2.0<-0.1)) 
    {fprintf(stderr,"triangulated surface is not closed. exit.\n"); return false;} 
  int edgeOrientation[numOfEdges][2];//edgeOrientation[i][0,1] stores from, to nodes. 
  for (int iTriangle=0; iTriangle<length_triangle_list; iTriangle++) {
    int from1, to1, from2, to2, from3, to3;
    bool found1, found2, found3;
    from1 = triangle_list[iTriangle][0];  to1 = triangle_list[iTriangle][1];  found1 = false;
    from2 = triangle_list[iTriangle][1];  to2 = triangle_list[iTriangle][2];  found2 = false;
    from3 = triangle_list[iTriangle][2];  to3 = triangle_list[iTriangle][0];  found3 = false;
    for (int j=0; j<iEdge; j++) {
      if ((edgeOrientation[j][0]==from1) && (edgeOrientation[j][1]==to1))
        {fprintf(stderr,"Not all triangles on SolidSurface point outward. exit.\n"); return false; }
      else if ((edgeOrientation[j][0]==to1) && (edgeOrientation[j][1]==from1))
        found1 = true;
      if ((edgeOrientation[j][0]==from2) && (edgeOrientation[j][1]==to2))
        {fprintf(stderr,"Not all triangles on SolidSurface point outward. exit.\n"); return false; }
      else if ((edgeOrientation[j][0]==to2) && (edgeOrientation[j][1]==from2))
        found2 = true;
      if ((edgeOrientation[j][0]==from3) && (edgeOrientation[j][1]==to3))
        {fprintf(stderr,"Not all triangles on SolidSurface point outward. exit.\n"); return false; }
      else if ((edgeOrientation[j][0]==to3) && (edgeOrientation[j][1]==from3))
        found3 = true;
    }
    if (!found1) {
      if (iEdge>=numOfEdges) 
        {fprintf(stderr,"triangulated surface is not closed. exit.\n"); return false;}
      edgeOrientation[iEdge][0] = from1;  edgeOrientation[iEdge][1] = to1; iEdge++;
    }
    if (!found2) {
      if (iEdge>=numOfEdges)
        {fprintf(stderr,"triangulated surface is not closed. exit.\n"); return false;}
      edgeOrientation[iEdge][0] = from2;  edgeOrientation[iEdge][1] = to2; iEdge++;
    }
    if (!found3) {
      if (iEdge>=numOfEdges)
        {fprintf(stderr,"triangulated surface is not closed. exit.\n"); return false;}
      edgeOrientation[iEdge][0] = from3;  edgeOrientation[iEdge][1] = to3; iEdge++;
    }
  }
  return true;
}

//---------------------------------------------------------------------------------

void DistEulerStructGhostFluid::constructInterface(int m, int n, int mn, double xmin, double xmax,
                                          double ymin, double ymax, double zmin, double zmax)
{
  fprintf(stderr,"DistEulerStructGhostFluid::constructInterface(m,n,mn,xmin,xmax...) called.\n");
  fprintf(stderr,"( m, n, mn ) = ( %d, %d, %d );  domain: ( %lf, %lf, %lf ) -> ( %lf, %lf, %lf )\n", 
                   m, n, mn, xmin, ymin, zmin, xmax, ymax, zmax);
  PhysBAM::GRID_3D<double> grid3d(m, n, mn, xmin, xmax, ymin, ymax, zmin, zmax);

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->constructInterface(grid3d);
}

//--------------------------------------------------------------------------------

void DistEulerStructGhostFluid::initializePhysBAMMPI(int numNodes,double depth)
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->initializePhysBAMMPI(numNodes,depth,com->cpuNum());
}

//--------------------------------------------------------------------------------

void DistEulerStructGhostFluid::computeLevelSet()
{
  fprintf(stderr,"DistEulerStructGhostFluid::computeLevelSet called.\n");

  // Initialize the Particles list
  PhysBAM::SOLIDS_PARTICLE<PhysBAM::VECTOR<double,3> >& physbam_solids_particle=*new PhysBAM::SOLIDS_PARTICLE<PhysBAM::VECTOR<double,3> >();
  physbam_solids_particle.Add_Particles(length_solids_particle_list);

  for (int i=0; i<length_solids_particle_list; i++) {
    physbam_solids_particle.X(i+1) = PhysBAM::VECTOR<double,3>(solids_particle_list[i][0], 
        solids_particle_list[i][1], solids_particle_list[i][2]);
  }

  std::cout<<physbam_solids_particle.X.Size()<<std::endl;

  // Initialize the Triangle list
  PhysBAM::LIST_ARRAY<PhysBAM::VECTOR<int,3> > & physbam_triangle_list=*new PhysBAM::LIST_ARRAY<PhysBAM::VECTOR<int,3> >();
  for (int i=0; i<length_triangle_list; i++){
    int nx, ny, nz;
    nx = triangle_list[i][0] + 1;  ny = triangle_list[i][1] + 1;  nz = triangle_list[i][2] + 1;
    physbam_triangle_list.Append(PhysBAM::VECTOR<int,3>(nx, ny, nz));
  }
  
  // Construct TRIANGLE_MESH triangle_mesh.
  PhysBAM::TRIANGLE_MESH& physbam_triangle_mesh=*new PhysBAM::TRIANGLE_MESH(physbam_solids_particle.number, physbam_triangle_list);

  // Construct TRIANGULATED_SURFACE.
  PhysBAM::TRIANGULATED_SURFACE<double>& physbam_triangulated_surface=*new PhysBAM::TRIANGULATED_SURFACE<double>(physbam_triangle_mesh, physbam_solids_particle);

//  PhysBAM_Interface->Write_Triangulated_Surface("triaSurf.tri",physbam_triangulated_surface);
  // Compute the Level Set
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->computeLevelSet(physbam_triangulated_surface);
}

//------------------------------------------------------------------------------
/*
bool DistEulerStructGhostFluid::insideOutside(double *position, const double xmin, const double xmax,
                     const double ymin, const double ymax, const double zmin, const double zmax)
{
  if (position[0] < xmin || position[0] > xmax || position[1] < ymin || position[1] > ymax ||
      position[2] < zmin || position[2] > zmax )  return false;
  return true;
}
*/
//-------------------------------------------------------------------------------

void DistEulerStructGhostFluid::computeGhostNodes()
{
  fprintf(stderr,"\n");
  double  tt = timer->getTime();

  fprintf(stderr,"DistEulerStructGhostFluid::computeGhostNodes called.\n" );
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->computeGhostNodes();

  fprintf(stderr,"\n");
  fprintf(stderr,"*** cpu no: %d: TIME for computeGhostNodes: %lf. ***\n", com->cpuNum(), timer->getTime()-tt);

/*
  if ((!PhysBAM_Interface)||(!solids_particle)||(!tetrahedralized_volume)||(!elem_list)||(!tetrahedron_mesh)) {
    fprintf(stderr,"ingredients lost!\n"); return;
  }
  PhysBAM_Interface->Compute_Ghost_Cells(*tetrahedralized_volume, *compressible_fluid_particle, bandwidth);
*/
}

//----------------------------------------------------------------------------------

void DistEulerStructGhostFluid::getMinAndMax(Vec3D* list, int &size, double *output) 
{
  for (int i=0; i<3; i++)  output[2*i] = output[2*i+1] = list[0][i];
  for (int iNode=0; iNode<size; iNode++) 
    for (int i=0; i<3; i++) {
      if (list[iNode][i]<output[2*i]) output[2*i] = list[iNode][i];
      if (list[iNode][i]>output[2*i+1]) output[2*i+1] = list[iNode][i];
    }
}

//-------------------------------------------------------------------------------------

void DistEulerStructGhostFluid::specifyBoundingBox(DistSVec<double,3> *X)
{
  if ((Xmax-Xmin)>1e-12 && (Ymax-Ymin)>1e-12 && (Zmax-Zmin)>1e-12) 
    fprintf(stderr,"WARNING:replacing existed bounding box.\n");
  double *mins = X->min(), *maxes = X->max();
  double percentX = (maxes[0] - mins[0])/100.0, percentY = (maxes[1] - mins[1])/100.0, percentZ = (maxes[2] - mins[2])/100.0;
  delete[] mins; delete[] maxes;
  if (!solids_particle_list) {fprintf(stderr,"Can't find Solid Surface."); return;}
  double minAndmax[6]; 
  getMinAndMax(solids_particle_list, length_solids_particle_list, minAndmax);
  fprintf(stderr,"Structure Surface: [%lf, %lf, %lf ] -> [%lf, %lf, %lf ] \n", minAndmax[0], minAndmax[2], minAndmax[4], minAndmax[1], minAndmax[3], minAndmax[5]);
  Xmin = minAndmax[0] - 10.0*percentX;  Xmax = minAndmax[1] + 10.0*percentX;
  Ymin = minAndmax[2] - 10.0*percentY;  Ymax = minAndmax[3] + 10.0*percentY;
  Zmin = minAndmax[4] - 10.0*percentZ;  Zmax = minAndmax[5] + 10.0*percentZ;
  fprintf(stderr,"BoundingBox: [%lf, %lf, %lf ] -> [%lf, %lf, %lf ] \n", Xmin,Ymin,Zmin,Xmax,Ymax,Zmax);
}

//--------------------------------------------------------------------------------------

double DistEulerStructGhostFluid::specifydx(Domain *domain, DistSVec<double,3>* X)
//given BoundingBox, find dx. 
{
  if ((Xmax-Xmin)<1e-12 && (Ymax-Ymin)<1e-12 && (Zmax-Zmin)<1e-12) {
    fprintf(stderr,"WARNING:BoundingBox undefined! return dx = 0.0.\n");
    return 0.0;
  }
  double minEdgeLength, aveEdgeLength, maxEdgeLength, dx;
  int numInsideEdges;
  domain->computeCharacteristicEdgeLength(*X, minEdgeLength, aveEdgeLength, maxEdgeLength, numInsideEdges,
                                          Xmin, Xmax, Ymin, Ymax, Zmin, Zmax);
  if (minEdgeLength < aveEdgeLength/4.0)  dx = aveEdgeLength/4.0;  else dx = minEdgeLength;
  
 
  fprintf(stderr,"minEdgeLength: %lf, aveEdgeLength: %lf, maxEdgeLength: %lf, numInsideEdges: %d, dx: %lf\n", minEdgeLength, aveEdgeLength, maxEdgeLength, numInsideEdges, dx);
  return dx;
}

//-------------------------------------------------------------------------------------

double DistEulerStructGhostFluid::specifyBandwidth()
{
  double subDomainBandwidth[numLocSub];
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) 
    subDomainBandwidth[iSub] = subESGF[iSub]->specifyBandwidth();

  double maxBandwidth = subDomainBandwidth[0];
  for (int iSub=0; iSub<numLocSub; iSub++)
    if (subDomainBandwidth[iSub]>maxBandwidth)  maxBandwidth = subDomainBandwidth[iSub];
  
  com->globalMax(1, &maxBandwidth);
  
  bandwidth = maxBandwidth*1.1;
  return bandwidth;
}

//--------------------------------------------------------------------------------------

bool DistEulerStructGhostFluid::updateStructureDynamics(double dt) //currently only for Forced Motion
{
  fprintf(stderr,"DistEulerStructGhostFluid::updateStructureDynamics(...) called.\n");
  if (vel.norm()>1.e-10) dsp += dt*vel; else return false;
  fprintf(stderr,"displacement within this time-step: (%f, %f, %f). dsp = (%f, %f, %f).\n", (dt*vel)[0], (dt*vel)[1], (dt*vel)[2], dsp[0], dsp[1], dsp[2]);
  return true;
}

//--------------------------------------------------------------------------------------

void DistEulerStructGhostFluid::recomputeLevelSet(Vec3D dsp)
{
  fprintf(stderr,"DistEulerStructGhostFluid::recomputeLevelSet(...) called.\n");
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->recomputeLevelSet(dsp);
}

//--------------------------------------------------------------------------------------

void DistEulerStructGhostFluid::updateRealNodeTag(DistSVec<double,3>& X, double xmin, double xmax, 
                                                  double ymin, double ymax, double zmin, double zmax)
{
  fprintf(stderr,"DistEulerStructGhostFluid::updateRealNodeTag() called.\n");
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->updateRealNodeTag(X(iSub), xmin, xmax, ymin, ymax, zmin, zmax);
} 








