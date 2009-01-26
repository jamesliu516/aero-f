#include <DistEulerStructGhostFluid.h>
#include <stdio.h>
//------------------------------------------------------------------------------------

DistEulerStructGhostFluid::DistEulerStructGhostFluid(Domain* domain, IoData& iod):numGlobNodes(domain->numNodes())
{
  numLocSub = domain->getNumLocSub();
  subDomain = domain->getSubDomain();

  philevel = new DistVec<double> (domain->getNodeDistInfo());
  gradPhilevel = new DistSVec<double,3>(domain->getNodeDistInfo());

  subESGF = new EulerStructGhostFluid *[numLocSub];
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub] = new EulerStructGhostFluid((*philevel)(iSub), (*gradPhilevel)(iSub),
                                               subDomain[iSub]->numNodes(), subDomain[iSub]->numElems());

  length_triangle_list = 0;
  length_solids_particle_list = 0;
  triangle_list = 0;
  solids_particle_list = 0;
  tag = 0;

  bandwidth = -1.0;
  totalForce[0] = totalForce[1] = totalForce[2] = 0.0;

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

  forceFile = fopen("/home/icmewang/Simulations/RiemanFSI/CPU16/results/force","w");
  pref = iod.ref.rv.pressure;

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

  if (philevel) delete philevel;
  if (gradPhilevel) delete gradPhilevel;

  if (solidsurface) delete[] solidsurface;

  fclose(forceFile);
}

//------------------------------------------------------------------------------------

void  DistEulerStructGhostFluid::prepareForCommunication()
{
  com->fprintf(stderr,"DistEulerStructGhostFluid::prepareForCommunication called. \n");
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

  com->fprintf(stderr,"solid surface: %d nodes, %d elements.\n", length_solids_particle_list, length_triangle_list);

  int thisNode;
  for (int iNode=0; iNode<length_solids_particle_list; iNode++)
    fscanf(topFile, "%d %lf %lf %lf", &thisNode, &(solids_particle_list[iNode][0]),
                                           &(solids_particle_list[iNode][1]), &(solids_particle_list[iNode][2]));
  if (thisNode!=length_solids_particle_list) {fprintf(stderr,"error in loading surface from file *!\n"); exit(1);}

  int nothing;
  for (int iElem=0; iElem<length_triangle_list; iElem++) {
    fscanf(topFile, "%d %d %d %d %d", &thisNode, &nothing, &(triangle_list[iElem][0]), &(triangle_list[iElem][1]),
                                         &(triangle_list[iElem][2]));
    triangle_list[iElem][0]--; triangle_list[iElem][1]--; triangle_list[iElem][2]--;
  }
  if (thisNode!=length_triangle_list) {fprintf(stderr,"error in loading surface from file **!\n", thisNode); exit(1);}
  fclose(topFile);

  //3. verify (1)triangulated surface is closed (2) normal's of all triangles point outward.
  com->fprintf(stderr,"Checking the solid surface...\n");
  if (checkTriangulatedSurface()) com->fprintf(stderr,"OK.\n");
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
  double mins[3], maxes[3];
  X->min(mins); X->max(maxes);
  double percentX = (maxes[0] - mins[0])/100.0, percentY = (maxes[1] - mins[1])/100.0, percentZ = (maxes[2] - mins[2])/100.0;
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

void DistEulerStructGhostFluid::initializePhysBAM(int m, int n, int mn, double xmin, double xmax,
                                          double ymin, double ymax, double zmin, double zmax)
{
  com->fprintf(stderr,"DistEulerStructGhostFluid::constructInterface(m,n,mn,xmin,xmax...) called.\n");
  com->fprintf(stderr,"( m, n, mn ) = ( %d, %d, %d );  domain: ( %lf, %lf, %lf ) -> ( %lf, %lf, %lf )\n",
                   m, n, mn, xmin, ymin, zmin, xmax, ymax, zmax);
  PhysBAM::GRID_3D<double> grid3d(m, n, mn, xmin, xmax, ymin, ymax, zmin, zmax);

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->initializePhysBAM(grid3d);
}

//--------------------------------------------------------------------------------

double DistEulerStructGhostFluid::specifyBandwidth()
{
  com->fprintf(stderr,"DistEulerStructGhostFluid::specifyBandwidth() called.\n");
  double subDomainBandwidth[numLocSub];
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomainBandwidth[iSub] = subDomain[iSub]->specifyBandwidth((*philevel)(iSub));

  double maxBandwidth = subDomainBandwidth[0];
  for (int iSub=0; iSub<numLocSub; iSub++)
    if (subDomainBandwidth[iSub]>maxBandwidth)  maxBandwidth = subDomainBandwidth[iSub];

  com->globalMax(1, &maxBandwidth);

  bandwidth = maxBandwidth*1.1;
  return bandwidth;
}

//--------------------------------------------------------------------------------

void DistEulerStructGhostFluid::initializePhysBAMMPI()
{
  com->fprintf(stderr,"DistEulerStructGhostFluid::initializePhysBAMMPI() called.\n");
  
  if (numGlobNodes<=0||bandwidth<=0.0) {fprintf(stderr,"ERROR: negative inputs. Aborting.\n");exit(-1);}
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) {
    int* locToGlobNodeMap = subDomain[iSub]->getNodeMap();
    subESGF[iSub]->initializePhysBAMMPI(locToGlobNodeMap,numGlobNodes,bandwidth);
  }
}

//------------------------------------------------------------------------------------

void DistEulerStructGhostFluid::computeLevelSet()
{
  com->fprintf(stderr,"DistEulerStructGhostFluid::computeLevelSet called.\n");

  // Initialize the Particles list
  PhysBAM::SOLIDS_PARTICLE<PhysBAM::VECTOR<double,3> >& physbam_solids_particle=*new PhysBAM::SOLIDS_PARTICLE<PhysBAM::VECTOR<double,3> >();
  physbam_solids_particle.Add_Particles(length_solids_particle_list);

  for (int i=0; i<length_solids_particle_list; i++) {
    physbam_solids_particle.X(i+1) = PhysBAM::VECTOR<double,3>(solids_particle_list[i][0],
        solids_particle_list[i][1], solids_particle_list[i][2]);
  }
//  std::cout<<physbam_solids_particle.X.Size()<<std::endl;

  // Initialize the Triangle list
  PhysBAM::LIST_ARRAY<PhysBAM::VECTOR<int,3> > & physbam_triangle_list=*new PhysBAM::LIST_ARRAY<PhysBAM::VECTOR<int,3> >();
  for (int i=0; i<length_triangle_list; i++){
    int nx, ny, nz;
    nx = triangle_list[i][0] + 1;  ny = triangle_list[i][1] + 1;  nz = triangle_list[i][2] + 1;
    physbam_triangle_list.Append(PhysBAM::VECTOR<int,3>(nx, ny, nz));
  }

  // Construct TRIANGLE_MESH triangle_mesh.
  PhysBAM::TRIANGLE_MESH& physbam_triangle_mesh=*new PhysBAM::TRIANGLE_MESH(physbam_solids_particle.number, physbam_triangle_list);
  physbam_triangle_mesh.Initialize_Adjacent_Elements();
  // Construct TRIANGULATED_SURFACE.
  PhysBAM::TRIANGULATED_SURFACE<double>& physbam_triangulated_surface=*new PhysBAM::TRIANGULATED_SURFACE<double>(physbam_triangle_mesh, physbam_solids_particle);
  physbam_triangulated_surface.Update_Triangle_List();

  // Compute the Level Set
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->computeLevelSet(physbam_triangulated_surface);
}

//-------------------------------------------------------------------------------

void DistEulerStructGhostFluid::getPhiFromModule(DistSVec<double, 3> &X, double xmin, double xmax,
                                                 double ymin, double ymax, double zmin, double zmax)
{
  com->fprintf(stderr, "DistEulerStructGhostFluid::getPhiFromModule called. \n");
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->getPhiFromModule(X(iSub), xmin, xmax, ymin, ymax, zmin, zmax);
}

//-------------------------------------------------------------------------------------

void DistEulerStructGhostFluid::getGradPhiFromModule(DistSVec<double, 3> &X, double xmin, double xmax,
                                                     double ymin, double ymax, double zmin, double zmax)
{
  com->fprintf(stderr,"DistEulerStructGhostFluid::getGradPhiFromModule called.\n");
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->getGradPhiFromModule(X(iSub), xmin, xmax, ymin, ymax, zmin, zmax);
}

//----------------------------------------------------------------------------------------

void DistEulerStructGhostFluid::getTetNearInterface(DistSVec<double,3> &X,
                                                    const double xmin, const double xmax, const double ymin,
                                                    const double ymax, const double zmin, const double zmax)
{
  com->fprintf(stderr, "DistEulerStructGhostFluid::getTetNearInterface called.\n");
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) {
    subDomain[iSub]->getMeshInBoundingBox(X(iSub), xmin, xmax, ymin, ymax, zmin, zmax, 
                                          subESGF[iSub]->nodeTag, subESGF[iSub]->numChosenNodes,
                                          subESGF[iSub]->tempNodeList, subESGF[iSub]->numChosenElems,
                                          subESGF[iSub]->tempElemList); 
    subESGF[iSub]->getTetNearInterface(X(iSub));
  }
}

//-----------------------------------------------------------------------------------

void DistEulerStructGhostFluid::clearTotalForce() 
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subESGF[iSub]->totalForce[0] = subESGF[iSub]->totalForce[1] = subESGF[iSub]->totalForce[2] = 0.0;
  totalForce[0] = totalForce[1] = totalForce[2] = 0.0;
}

//-----------------------------------------------------------------------------------

Vec3D DistEulerStructGhostFluid::getTotalForce()
{
  double subDomainTotalForce[numLocSub][3];
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) {
    subDomainTotalForce[iSub][0] = subESGF[iSub]->totalForce[0];
    subDomainTotalForce[iSub][1] = subESGF[iSub]->totalForce[1];
    subDomainTotalForce[iSub][2] = subESGF[iSub]->totalForce[2];
  }
  for (int iSub=0; iSub<numLocSub; iSub++) {
    totalForce[0] += subDomainTotalForce[iSub][0];
    totalForce[1] += subDomainTotalForce[iSub][1];
    totalForce[2] += subDomainTotalForce[iSub][2];
  }
  com->globalSum(3, totalForce);

  com->fprintf(stderr,"Total Force on Structure Surface = [%f, %f, %f]\n", totalForce[0]*pref, totalForce[1]*pref, totalForce[2]*pref);
  com->fprintf(forceFile, "%lf %lf %lf\n", totalForce[0]*pref, totalForce[1]*pref, totalForce[2]*pref);
}



























