#include <EulerStructGhostFluid.h>

//---------------------------------------------------------------------------------------
/*template<int dim>
void EulerStructGhostFluid::mirroring(SVec<double,dim>& V)
{
  for (int iNode=0; iNode<numGNode; iNode++){
    int* elemNodeNum = subDomain->getElemNodeNum(tetList[iNode]);

//Special case 1.
    if (elemNodeNum[0] == gNodeList[iNode]){
      double alpha, beta, gamma;
      alpha = localCoord[iNode][0];  beta = localCoord[iNode][1];  gamma = localCoord[iNode][2];

      //now compute rho and p
      double rho1, rho2, rho3, rho;
      rho1 = V[elemNodeNum[1]][0];  rho2 = V[elemNodeNum[2]][0];  rho3 = V[elemNodeNum[3]][0];
      rho = (alpha*rho1 + beta*rho2 + gamma*rho3) / (alpha + beta + gamma);            
      double p1, p2, p3, p;
      p1 = V[elemNodeNum[1]][4];  p2 = V[elemNodeNum[2]][4];  p3 = V[elemNodeNum[3]][4];
      p = (alpha*p1 + beta*p2 + gamma*p3) / (alpha + beta + gamma);

      //now compute u, v, w.
      double c = 1.0 - alpha - beta - gamma;
      double xi1, xi2, xi3;
      xi1 = alpha*V[elemNodeNum[1]][1] + beta*V[elemNodeNum[2]][1] + gamma*V[elemNodeNum[3]][1];
      xi2 = alpha*V[elemNodeNum[1]][2] + beta*V[elemNodeNum[2]][2] + gamma*V[elemNodeNum[3]][2];
      xi3 = alpha*V[elemNodeNum[1]][3] + beta*V[elemNodeNum[2]][3] + gamma*V[elemNodeNum[3]][3];
      
      double n1, n2, n3;
      n1 = n_phi[iNode][0];  n2 = n_phi[iNode][1];  n3 = n_phi[iNode][2];
      
      double a11, a12, a13, a21, a22, a23, a31, a32, a33;
      a11 = 1.0-c+2.0*c*n1*n1;  a12 = 2.0*c*n1*n2;        a13 = 2.0*c*n1*n3;
      a21 = 2.0*c*n2*n1;        a22 = 1.0-c+2.0*c*n2*n2;  a23 = 2.0*c*n2*n3;
      a31 = 2.0*c*n3*n1;        a32 = 2.0*c*n3*n2;        a33 = 1.0-c+2.0*c*n3*n3; 

      double u, v, w;
      solve3by3LinearSystem(a11, a12, a13, a21, a22, a23, a31, a32, a33, 
                            u, v, w, xi1, xi2, xi3);
      double u_n = u*n1 + v*n2 + w*n3;
      u = u - 2.0*u_n*n1;  v = v - 2.0*u_n*n2;  w = w -2.0*u_n*n3;
    
      V[iNode][0] = rho;
      V[iNode][1] = u;
      V[iNode][2] = v;
      V[iNode][3] = w;
      V[iNode][4] = p;

      continue;
    }



//Special case 2.
    if (elemNodeNum[1] == gNodeList[iNode]){
      double alpha, beta, gamma;
      alpha = localCoord[iNode][0];  beta = localCoord[iNode][1];  gamma = localCoord[iNode][2];

      //now compute rho and p.
      double rho0, rho2, rho3, rho;
      rho0 = V[elemNodeNum[0]][0];  rho2 = V[elemNodeNum[2]][0];  rho3 = V[elemNodeNum[3]][0];
      rho = (rho0 - alpha*rho0 + beta*(rho2 - rho0) + gamma*(rho3 - rho0)) / (1.0 - alpha);
      double p0, p2, p3, p;
      p0 = V[elemNodeNum[0]][4];  p2 = V[elemNodeNum[2]][4];  p3 = V[elemNodeNum[3]][4];
      p = (p0 - alpha*p0 + beta*(p2 - p0) + gamma*(p3 - p0)) / (1.0 - alpha);
     
      //now compute u, v, w.
      double c = alpha;
      double xi1, xi2, xi3;
      xi1 = (1.0-alpha)*V[elemNodeNum[0]][1] 
            + beta*(V[elemNodeNum[2]][1] - V[elemNodeNum[0]][1])
            + gamma*(V[elemNodeNum[3]][1] - V[elemNodeNum[0]][1]);
      xi2 = (1.0-alpha)*V[elemNodeNum[0]][2] 
            + beta*(V[elemNodeNum[2]][2] - V[elemNodeNum[0]][2])
            + gamma*(V[elemNodeNum[3]][2] - V[elemNodeNum[0]][2]);
      xi3 = (1.0-alpha)*V[elemNodeNum[0]][3]
            + beta*(V[elemNodeNum[2]][3] - V[elemNodeNum[0]][3])
            + gamma*(V[elemNodeNum[3]][3] - V[elemNodeNum[0]][3]);

      double n1, n2, n3;
      n1 = n_phi[iNode][0];  n2 = n_phi[iNode][1];  n3 = n_phi[iNode][2];

      double a11, a12, a13, a21, a22, a23, a31, a32, a33;
      a11 = 1.0-c+2.0*c*n1*n1;  a12 = 2.0*c*n1*n2;        a13 = 2.0*c*n1*n3;
      a21 = 2.0*c*n2*n1;        a22 = 1.0-c+2.0*c*n2*n2;  a23 = 2.0*c*n2*n3;
      a31 = 2.0*c*n3*n1;        a32 = 2.0*c*n3*n2;        a33 = 1.0-c+2.0*c*n3*n3;

      double u, v, w;
      solve3by3LinearSystem(a11, a12, a13, a21, a22, a23, a31, a32, a33,
                            u, v, w, xi1, xi2, xi3);
      double u_n = u*n1 + v*n2 + w*n3;
      u = u - 2.0*u_n*n1;  v = v - 2.0*u_n*n2;  w = w -2.0*u_n*n3;

      V[iNode][0] = rho;
      V[iNode][1] = u;
      V[iNode][2] = v;
      V[iNode][3] = w;
      V[iNode][4] = p;

      continue;
    }



//Special case 2.
    if (elemNodeNum[2] == gNodeList[iNode]){
      double alpha, beta, gamma;
      alpha = localCoord[iNode][0];  beta = localCoord[iNode][1];  gamma = localCoord[iNode][2];

      //now compute rho and p.
      double rho0, rho1, rho3, rho;
      rho0 = V[elemNodeNum[0]][0];  rho1 = V[elemNodeNum[1]][0];  rho3 = V[elemNodeNum[3]][0];
      rho = (rho0 - beta*rho0 + alpha*(rho1 - rho0) + gamma*(rho3 - rho0)) / (1.0 - beta);
      double p0, p1, p3, p;
      p0 = V[elemNodeNum[0]][4];  p1 = V[elemNodeNum[2]][4];  p3 = V[elemNodeNum[3]][4];
      p = (p0 - beta*p0 + alpha*(p1 - p0) + gamma*(p3 - p0)) / (1.0 - beta);

      //now compute u, v, w.
      double c = beta;
      double xi1, xi2, xi3;
      xi1 = (1.0-beta)*V[elemNodeNum[0]][1]
            + alpha*(V[elemNodeNum[1]][1] - V[elemNodeNum[0]][1])
            + gamma*(V[elemNodeNum[3]][1] - V[elemNodeNum[0]][1]);
      xi2 = (1.0-beta)*V[elemNodeNum[0]][2]
            + alpha*(V[elemNodeNum[1]][2] - V[elemNodeNum[0]][2])
            + gamma*(V[elemNodeNum[3]][2] - V[elemNodeNum[0]][2]);
      xi3 = (1.0-beta)*V[elemNodeNum[0]][3]
            + alpha*(V[elemNodeNum[1]][3] - V[elemNodeNum[0]][3])
            + gamma*(V[elemNodeNum[3]][3] - V[elemNodeNum[0]][3]);

      double n1, n2, n3;
      n1 = n_phi[iNode][0];  n2 = n_phi[iNode][1];  n3 = n_phi[iNode][2];

      double a11, a12, a13, a21, a22, a23, a31, a32, a33;
      a11 = 1.0-c+2.0*c*n1*n1;  a12 = 2.0*c*n1*n2;        a13 = 2.0*c*n1*n3;
      a21 = 2.0*c*n2*n1;        a22 = 1.0-c+2.0*c*n2*n2;  a23 = 2.0*c*n2*n3;
      a31 = 2.0*c*n3*n1;        a32 = 2.0*c*n3*n2;        a33 = 1.0-c+2.0*c*n3*n3;

      double u, v, w;
      solve3by3LinearSystem(a11, a12, a13, a21, a22, a23, a31, a32, a33,
                            u, v, w, xi1, xi2, xi3);
      double u_n = u*n1 + v*n2 + w*n3;
      u = u - 2.0*u_n*n1;  v = v - 2.0*u_n*n2;  w = w -2.0*u_n*n3;

      V[iNode][0] = rho;
      V[iNode][1] = u;
      V[iNode][2] = v;
      V[iNode][3] = w;
      V[iNode][4] = p;

      continue;
    }


//Special case 2.
    if (elemNodeNum[3] == gNodeList[iNode]){
      double alpha, beta, gamma;
      alpha = localCoord[iNode][0];  beta = localCoord[iNode][1];  gamma = localCoord[iNode][2];

      //now compute rho and p.
      double rho0, rho1, rho2, rho;
      rho0 = V[elemNodeNum[0]][0];  rho1 = V[elemNodeNum[1]][0];  rho2 = V[elemNodeNum[2]][0];
      rho = (rho0 - gamma*rho0 + alpha*(rho1 - rho0) + beta*(rho2 - rho0)) / (1.0 - gamma);
      double p0, p1, p2, p;
      p0 = V[elemNodeNum[0]][4];  p1 = V[elemNodeNum[2]][4];  p2 = V[elemNodeNum[2]][4];
      p = (p0 - gamma*p0 + alpha*(p1 - p0) + beta*(p2 - p0)) / (1.0 - gamma);

      //now compute u, v, w.
      double c = gamma;
      double xi1, xi2, xi3;
      xi1 = (1.0-gamma)*V[elemNodeNum[0]][1]
            + alpha*(V[elemNodeNum[1]][1] - V[elemNodeNum[0]][1])
            + beta*(V[elemNodeNum[2]][1] - V[elemNodeNum[0]][1]);
      xi2 = (1.0-gamma)*V[elemNodeNum[0]][2]
            + alpha*(V[elemNodeNum[1]][2] - V[elemNodeNum[0]][2])
            + beta*(V[elemNodeNum[2]][2] - V[elemNodeNum[0]][2]);
      xi3 = (1.0-gamma)*V[elemNodeNum[0]][3]
            + alpha*(V[elemNodeNum[1]][3] - V[elemNodeNum[0]][3])
            + beta*(V[elemNodeNum[2]][3] - V[elemNodeNum[0]][3]);

      double n1, n2, n3;
      n1 = n_phi[iNode][0];  n2 = n_phi[iNode][1];  n3 = n_phi[iNode][2];

      double a11, a12, a13, a21, a22, a23, a31, a32, a33;
      a11 = 1.0-c+2.0*c*n1*n1;  a12 = 2.0*c*n1*n2;        a13 = 2.0*c*n1*n3;
      a21 = 2.0*c*n2*n1;        a22 = 1.0-c+2.0*c*n2*n2;  a23 = 2.0*c*n2*n3;
      a31 = 2.0*c*n3*n1;        a32 = 2.0*c*n3*n2;        a33 = 1.0-c+2.0*c*n3*n3;

      double u, v, w;
      solve3by3LinearSystem(a11, a12, a13, a21, a22, a23, a31, a32, a33,
                            u, v, w, xi1, xi2, xi3);
      double u_n = u*n1 + v*n2 + w*n3;
      u = u - 2.0*u_n*n1;  v = v - 2.0*u_n*n2;  w = w -2.0*u_n*n3;

      V[iNode][0] = rho;
      V[iNode][1] = u;
      V[iNode][2] = v;
      V[iNode][3] = w;
      V[iNode][4] = p;

      continue;
    }



//General case.
    double* Vtemp = V[elemNodeNum[0]] + localCoord[iNode][0]*(V[elemNodeNum[1]]-V[elemNodeNum[0]])
                                     + localCoord[iNode][1]*(V[elemNodeNum[2]]-V[elemNodeNum[0]])
                                     + localCoord[iNode][2]*(V[elemNodeNum[3]]-V[elemNodeNum[0]]);
  
    //now we update V[iNode].
    double rho, u, v, w, p;
    rho = Vtemp[0];  u = Vtemp[1];  v = Vtemp[2];  w = Vtemp[3];  p = Vtemp[4];
    double u_n = u*n_phi[iNode][0] + v*n_phi[iNode][1] + w*n_phi[iNode][2];
    u = u - 2.0*u_n*n_phi[iNode][0];
    v = v - 2.0*u_n*n_phi[iNode][1];
    w = w - 2.0*u_n*n_phi[iNode][2];
    
    //now need to change V[iNode] = {rho, u, v, w, p};
    V[iNode][0] = rho;
    V[iNode][1] = u;
    V[iNode][2] = v;
    V[iNode][3] = w;
    V[iNode][4] = p;


  }
}

*/

//-------------------------------------------------------------------------------------------------
template<int dim>
void EulerStructGhostFluid::getTetNearInterface(SVec<double,3> &subX, 
                                                const double xmin, const double xmax, const double ymin,
                                                const double ymax, const double zmin, const double zmax,
                                                SVec<double,dim> &subU, const double bandwidth, 
                                                bool skipEulerianGhostInterpolation = false)
{
  //3. build a temporary element list.
  int numElems = subDomain->numElems();
  if (tempElemList) delete[] tempElemList;  tempElemList = new int[numElems];
  numChosenElems = 0;

  for (int i=0; i<numElems; i++) {
    int* elemNodeNum = subDomain->getElemNodeNum(i);
    // check if any node is (1)inside bounding box & (2) within the "band". 
    if ((this->insideOutside(subX[elemNodeNum[0]], xmin, xmax, ymin, ymax, zmin, zmax) ||
        this->insideOutside(subX[elemNodeNum[1]], xmin, xmax, ymin, ymax, zmin, zmax)  ||
        this->insideOutside(subX[elemNodeNum[2]], xmin, xmax, ymin, ymax, zmin, zmax)  ||
        this->insideOutside(subX[elemNodeNum[3]], xmin, xmax, ymin, ymax, zmin, zmax)   ) 
//       ((philevel[elemNodeNum[0]] < bandwidth  &&  philevel[elemNodeNum[0]] > -bandwidth) ||
//        (philevel[elemNodeNum[1]] < bandwidth  &&  philevel[elemNodeNum[1]] > -bandwidth) ||
//        (philevel[elemNodeNum[2]] < bandwidth  &&  philevel[elemNodeNum[2]] > -bandwidth) ||
//        (philevel[elemNodeNum[3]] < bandwidth  &&  philevel[elemNodeNum[3]] > -bandwidth))) 
       )
/*
       ((philevel[elemNodeNum[0]] > -bandwidth) ||
        (philevel[elemNodeNum[1]] > -bandwidth) ||
        (philevel[elemNodeNum[2]] > -bandwidth) ||
        (philevel[elemNodeNum[3]] > -bandwidth)))
*/
    {
      tempElemList[numChosenElems] = i;
      numChosenElems++;
    }
  }
    //4. construct a tag on nodes. build a temporary node list.
  int numNodes = subDomain->numNodes();
  if (tempNodeList) delete[] tempNodeList;  tempNodeList = new int[numNodes];
  if (nodeTag) delete[] nodeTag;  nodeTag = new int[numNodes];
  numChosenNodes = 0;
  for (int i=0; i<numNodes; i++) nodeTag[i] = -1;

  for (int iElem=0; iElem<numChosenElems; iElem++) {
    int* elemNodeNum = subDomain->getElemNodeNum(tempElemList[iElem]);
    for (int j=0; j<4; j++)
      if (nodeTag[elemNodeNum[j]] < 0) {
        nodeTag[elemNodeNum[j]] = numChosenNodes;
        tempNodeList[numChosenNodes] = elemNodeNum[j];
        numChosenNodes++;
      }
  }
  //5.construct Compressible_Fluid_Particle_List.
  if (compressible_fluid_particle) delete compressible_fluid_particle;
  compressible_fluid_particle = new PhysBAM::COMPRESSIBLE_FLUID_PARTICLE<PhysBAM::VECTOR<double,3> > ;
  compressible_fluid_particle->Add_Particles(numChosenNodes);
  for (int iParticle=0; iParticle<numChosenNodes; iParticle++) {
    if (nodeTag[tempNodeList[iParticle]] != iParticle) {
      fprintf(stderr,"error in constructing COMPRESSIBLE_FLUID_PARTICLE!\n");
      exit(1);
    }
    double* position = subX[tempNodeList[iParticle]];
    compressible_fluid_particle->X(iParticle+1) = PhysBAM::VECTOR<double,3>(position[0], position[1], position
[2]);

    compressible_fluid_particle->rho(iParticle+1) = subU[tempNodeList[iParticle]][0];
    double v1 = subU[tempNodeList[iParticle]][1];
    double v2 = subU[tempNodeList[iParticle]][2];
    double v3 = subU[tempNodeList[iParticle]][3];
    compressible_fluid_particle->V(iParticle+1) = PhysBAM::VECTOR<double,3>(v1, v2, v3);
    compressible_fluid_particle->E(iParticle+1) = subU[tempNodeList[iParticle]][4];
//      if (iParticle<=10) fprintf(stderr,"Node: %d, values: %lf, (%lf, %lf, %lf), %lf.\n", iParticle+1, compressible_fluid_particle->rho(iParticle+1), v1, v2, v3, compressible_fluid_particle->E(iParticle+1));

  }
  //6.construct SOLIDS_PARTICLE.
  if (solids_particle) delete solids_particle;
  solids_particle = new PhysBAM::SOLIDS_PARTICLE<PhysBAM::VECTOR<double,3> >;
  solids_particle->Add_Particles(numChosenNodes);
  for (int iParticle=0; iParticle<numChosenNodes; iParticle++) {
    double* position = subX[tempNodeList[iParticle]];
    solids_particle->X(iParticle+1) = PhysBAM::VECTOR<double,3>(position[0], position[1], position[2]);
  }

  //7.construct LIST_ARRAY<VECTOR<int,4> >. (tetrahedron_list)
  if (elem_list) delete elem_list;
  elem_list = new PhysBAM::LIST_ARRAY<PhysBAM::VECTOR<int,4> >;
  for (int iElem=0; iElem<numChosenElems; iElem++) {
    int* elemNodeNum = subDomain->getElemNodeNum(tempElemList[iElem]);
    elem_list->Append(PhysBAM::VECTOR<int,4>(nodeTag[elemNodeNum[0]]+1, nodeTag[elemNodeNum[1]]+1, nodeTag[elemNodeNum[2]]+1, nodeTag[elemNodeNum[3]]+1));
  }

  //8.construct TETRAHEDRON_MESH and then TETRAHEDRALIZED_VOLUME.
  if (tetrahedron_mesh) delete tetrahedron_mesh;
  tetrahedron_mesh = new PhysBAM::TETRAHEDRON_MESH (numChosenNodes, *elem_list);
  if (tetrahedralized_volume) delete tetrahedralized_volume;
  tetrahedralized_volume = new PhysBAM::TETRAHEDRALIZED_VOLUME<double> (*tetrahedron_mesh, *solids_particle);
//  PhysBAM_Interface->Initialize_Acceleration_Structures(*tetrahedralized_volume,bandwidth);
  PhysBAM_Interface->Set_Compute_ALE_Ghost_Cells_Directly(skipEulerianGhostInterpolation);

  //For debugging only.
//  fprintf(stderr,"total elems = %d,  total nodes = %d.\n", subDomain[iSub]->numElems(), subDomain[iSub]->numNodes());
  fprintf(stderr,"numChosenElems = %d,  numChosenNodes = %d.\n", numChosenElems, numChosenNodes);

}

//------------------------------------------------------------------------------------

template<int dim>
void EulerStructGhostFluid::updateGhostNodes(SVec<double,dim> &subU, const double bandwidth_to_update,
                                             Vec3D vel)
{
  fprintf(stderr,"input velocity is (%f, %f, %f).\n", vel[0], vel[1], vel[2]);
  for (int iParticle=0; iParticle<numChosenNodes; iParticle++) 
    if (philevel[tempNodeList[iParticle]] >= -bandwidth_to_update ) {
      PhysBAM::VECTOR<double,3> velocity = compressible_fluid_particle->V(iParticle+1);
      
      subU.v[tempNodeList[iParticle]][0] = compressible_fluid_particle->rho(iParticle+1);
      subU.v[tempNodeList[iParticle]][1] = velocity(1) + 2.0*vel[0]*subU.v[tempNodeList[iParticle]][0];
      subU.v[tempNodeList[iParticle]][2] = velocity(2) + 2.0*vel[1]*subU.v[tempNodeList[iParticle]][0];
      subU.v[tempNodeList[iParticle]][3] = velocity(3) + 2.0*vel[2]*subU.v[tempNodeList[iParticle]][0];
      subU.v[tempNodeList[iParticle]][4] = compressible_fluid_particle->E(iParticle+1);

      if (subU[iParticle][0] < 0.f) {
        fprintf(stderr, "Error: got negative pressure on ghost nodes!\n"); 
        exit(0);
      }
      if (subU[iParticle][4] < 0.f) {
       fprintf(stderr, "Error: got negative energy (total energy) on ghost nodes!\n");
       exit(0);
      }
      if (0.4*(subU[iParticle][4] - 0.5*(subU[iParticle][1]*subU[iParticle][1] + 
                                     subU[iParticle][2]*subU[iParticle][2] +
                                     subU[iParticle][3]*subU[iParticle][3])/(subU[iParticle][0])) < 0.f){
        fprintf(stderr, "Error: got negative pressure on ghost nodes!\n");
        exit(0);
      }
    } 
}

//------------------------------------------------------------------------------------

template<int dim>
void EulerStructGhostFluid::updateCFP(SVec<double,dim> &subU, int cpuNum, int cpuSize)
{
  int count = 0;
  if (!compressible_fluid_particle) {fprintf(stderr,"compressible_fluid_particle doesn't exist!\n"); return;}
  for (int iParticle=0; iParticle<numChosenNodes; iParticle++) {
    count++;
    compressible_fluid_particle->rho(iParticle+1) = subU[tempNodeList[iParticle]][0];
    double v1 = subU[tempNodeList[iParticle]][1];
    double v2 = subU[tempNodeList[iParticle]][2];
    double v3 = subU[tempNodeList[iParticle]][3];
    compressible_fluid_particle->V(iParticle+1) = PhysBAM::VECTOR<double,3>(v1, v2, v3);
    compressible_fluid_particle->E(iParticle+1) = subU[tempNodeList[iParticle]][4];
  }
  fprintf(stderr,"cpu no: %d, subdomain no: %d, numParticles: %d, numCPU: %d.\n", cpuNum, subDomain->getLocSubNum(), count, cpuSize);

}

//--------------------------------------------------------------------------------------

template<int dim>
void EulerStructGhostFluid::computeTotalForce(Vec3D &forceToWrite, int (*triangles)[3], int lenTriangles, Vec3D* particles , double lenParticles, SVec<double,3>& subX, SVec<double,dim>& subU)
{
  fprintf(stderr,"EulerStructGhostFluid::computeTotalForce(...) called. num of Triangles: %d.\n",lenTriangles);
  forceToWrite = 0.0;
  if (chainMail) delete[] chainMail;  chainMail = new int[lenTriangles];
  int count = 0;

  for (int iTriangle=0; iTriangle<lenTriangles; iTriangle++) {
//    fprintf(stderr,"iTriangle = %d.\n", iTriangle);
    Vec3D particle0(particles[triangles[iTriangle][0]][0], particles[triangles[iTriangle][0]][1],
                    particles[triangles[iTriangle][0]][2]);
    Vec3D particle1(particles[triangles[iTriangle][1]][0], particles[triangles[iTriangle][1]][1],
                    particles[triangles[iTriangle][1]][2]);
    Vec3D particle2(particles[triangles[iTriangle][2]][0], particles[triangles[iTriangle][2]][1],
                    particles[triangles[iTriangle][2]][2]);

    //1. Find the centroid.
    Vec3D centroid = (particle0 + particle1 + particle2)/3.0;
    //2. Find the area.
    Vec3D normal = (particle1 - particle0)^(particle2 - particle0);
    Vec3D normal1 = normal;
    double area = 0.5*normal.norm();
    //3. Find the unit normal.
    normal = -1.0/normal.norm()*normal;
    //4. find the surranding tet. 
    int myTet = -1;
    for (int iElem=0; iElem<numChosenElems; iElem++) {
      int thisElem = tempElemList[iElem];
      int* thisElemNodes = subDomain->getElemNodeNum(thisElem);
//      if (philevel[thisElemNodes[0]]<-1e-8 && philevel[thisElemNodes[1]]<-1e-8 &&
//          philevel[thisElemNodes[2]]<-1e-8 && philevel[thisElemNodes[3]]<-1e-8) continue;

      if (subDomain->isINodeinITet(centroid, thisElem, subX))  {myTet = thisElem; break;}
    }

    if (myTet<0) {
      fprintf(stderr,"myTet = %d. Failed in finding the surrounding tet for Triange #%d. exit.\n", myTet, iTriangle);
      exit(-1);
    }
    chainMail[iTriangle] = myTet;

    // for debugging. 
/*    int* myTetNodes = subDomain->getElemNodeNum(myTet);
    if (philevel[myTetNodes[0]]>0.0 && philevel[myTetNodes[1]]>0.0 && 
        philevel[myTetNodes[2]]>0.0 && philevel[myTetNodes[3]]>0.0)
    {
      count++;
      fprintf(stderr, "Triangle #%d, A[%f, %f, %f], B[%f, %f, %f], C[%f, %f, %f].\n", iTriangle, particle0[0], particle0[1], particle0[2], particle1[0], particle1[1], particle1[2], particle2[0], particle2[1], particle2[2]);
      fprintf(stderr,"centroid: [%f, %f, %f,].\n", centroid[0], centroid[1], centroid[2]);
      fprintf(stderr,"tet nodes T1[%f, %f, %f], T2[%f, %f, %f], T3[%f, %f, %f], T4[%f, %f, %f].\n", subX[myTetNodes[0]][0], subX[myTetNodes[0]][1], subX[myTetNodes[0]][2], subX[myTetNodes[1]][0], subX[myTetNodes[1]][1], subX[myTetNodes[1]][2], subX[myTetNodes[2]][0], subX[myTetNodes[2]][1], subX[myTetNodes[2]][2], subX[myTetNodes[3]][0], subX[myTetNodes[3]][1], subX[myTetNodes[3]][2]);
      fprintf(stderr,"levelset ls(T1) = %f, ls(T2) = %f, ls(T3) = %f, ls(T4) = %f.\n", philevel[myTetNodes[0]], philevel[myTetNodes[1]], philevel[myTetNodes[2]], philevel[myTetNodes[3]]);
      exit(-1);
    }
    if (philevel[myTetNodes[0]]<0.0 && philevel[myTetNodes[1]]<0.0 &&
        philevel[myTetNodes[2]]<0.0 && philevel[myTetNodes[3]]<0.0)
      count++;
*/

    //5. get pressure at centroid by constant extrapolation.
    double pressure = -1.0;

    Vec3D normalLS = getGradPhiAtPosition(centroid);
    int* myNodes = subDomain->getElemNodeNum(myTet);
    double pCandidates[4];

    for (int iNode=0; iNode<4; iNode++)
      pCandidates[iNode] = 0.4*(subU[myNodes[iNode]][4] - 0.5*(subU[myNodes[iNode]][1]*
                           subU[myNodes[iNode]][1]+subU[myNodes[iNode]][2]*subU[myNodes[iNode]][2] +
                           subU[myNodes[iNode]][3]*subU[myNodes[iNode]][3])/(subU[myNodes[iNode]][0]));

    pressure = subDomain->scalarNormalExtrap(pCandidates, centroid, normalLS, myTet, subX, 1);
    pressure = pressure*1.0e5/7.936508;

    if (centroid[1]>0) pressure = 25.0; else pressure = 15.0;

// (so-called) extrapolation from closest node
//    int* myNodes = subDomain->getElemNodeNum(myTet);
//    for (int iNode=0; iNode<4; iNode++)
//      if (philevel[myNodes[iNode]]>0) {    //real node?
//        pressure = 0.4*(subU[myNodes[iNode]][4] - 0.5*(subU[myNodes[iNode]][1]*subU[myNodes[iNode]][1] +
//                        subU[myNodes[iNode]][2]*subU[myNodes[iNode]][2] +
//                        subU[myNodes[iNode]][3]*subU[myNodes[iNode]][3])/(subU[myNodes[iNode]][0]));
//        pressure = pressure*1.0e5/7.936508;
//      //fprintf(stderr,"iTriangle = %d, pressure = %e.\n",iTriangle, pressure);
//        break;
//      }

    if (pressure<=0.0) {
      fprintf(stderr,"negative pressure at Triangle #%d. exit.\n", iTriangle);
      exit(-1);
    }

    //6. update total force.
    Vec3D locForce = pressure*area*normal;
    forceToWrite += locForce;
  }
  fprintf(stderr,"number of centroids covered by real/ghost tets: %d.\n", count);
}

//--------------------------------------------------------------------------------------------------

template<int dim>
void EulerStructGhostFluid::populateNewFluid(SVec<double,3> &X, SVec<double,dim> &U)
{
  for (int iNode=0; iNode<subDomain->numNodes(); iNode++)
    if ((!realNodeTag0[iNode])&&realNodeTag[iNode]) {
      // arrived at a new real node
      fprintf(stderr,"found new real node: #%d [%f, %f, %f]. U[%d] = [%f, %f, %f, %f, %f].\n", 
              iNode, X[iNode][0], X[iNode][1], X[iNode][2], iNode, U[iNode][0], U[iNode][1],
              U[iNode][2], U[iNode][3], U[iNode][4]);
      Vec3D coords(X[iNode][0], X[iNode][1], X[iNode][2]);
      // get neighbour elements.
      int numNbElems;
      int* nbElems = subDomain->getNeiElemOfNode(iNode,1,numNbElems);
      fprintf(stderr,"number of neighbor elems: %d.\n", numNbElems);
      // get gradient of levelset at node 
      Vec3D normalLS = getGradPhiAtPosition(coords);
      normalLS = 1.0/normalLS.norm()*normalLS;

      // extrapolate along the gradient of levelset
      double tempU[dim];    
      int howmany = 0;
      bool found = false;
      for (int iElem=0; iElem<numNbElems; iElem++) {
        int* myNodes = subDomain->getElemNodeNum(nbElems[iElem]);
        int otherNodes[3];
        int count = 0; 

        for (int i=0; i<4; i++) 
          if (myNodes[i]!=iNode) {otherNodes[count]=myNodes[i]; count++;}
        if (count!=3) {fprintf(stderr,"error! myNodes = [%d, %d, %d, %d], iNode = %d, count = %d.\n", myNodes[0], myNodes[1], myNodes[2], myNodes[3], iNode, count); exit(-1);}

        Vec3D A(X[otherNodes[0]][0], X[otherNodes[0]][1], X[otherNodes[0]][2]);
        Vec3D B(X[otherNodes[1]][0], X[otherNodes[1]][1], X[otherNodes[1]][2]);
        Vec3D C(X[otherNodes[2]][0], X[otherNodes[2]][1], X[otherNodes[2]][2]);
//        fprintf(stderr,"A[%f, %f, %f], B[%f, %f, %f], C[%f, %f, %f]", A[0], A[1], A[2], 
//                       B[0], B[1], B[2], C[0], C[1], C[2]);
        Vec3D faceNormal = (B-A)^(C-A);  faceNormal = 1.0/faceNormal.norm()*faceNormal;
        double r = ((A-coords)*faceNormal)/(faceNormal*normalLS);
        if (r<0) continue; // opposite direction.

        Vec3D P = coords + r*normalLS;

        //get barycentric coords.
        double rA = 0.5*((B-P)^(C-P)).norm();
        double rB = 0.5*((C-P)^(A-P)).norm();
        double rC = 0.5*((A-P)^(B-P)).norm();
        double rU = 0.5*((B-A)^(C-A)).norm();
        if (rA<0.0) rA=-rA;  if (rB<0.0) rB=-rB;  if (rC<0.0) rC=-rC;  if (rU<0.0) rU=-rU;
        rA /= rU;  rB /= rU;  rC /= rU;

//        fprintf(stderr,"rA = %e, rB = %e, rC = %e. sum-1 = %e\n", rA, rB, rC, rA+rB+rC-1.0);
        if ((rA+rB+rC-1.0)>1e-10 || (rA+rB+rC-1.0)<-1e-10) continue;
        found = true; howmany++;

        // update tempU
        for (int i=0; i<dim; i++)
          tempU[i] = U[otherNodes[0]][i]*rA + U[otherNodes[1]][i]*rB + U[otherNodes[2]][i]*rC;
        // check if data come from real nodes.
        if (philevel[otherNodes[0]]<0 || philevel[otherNodes[1]]<0 || philevel[otherNodes[2]]<0)
          fprintf(stderr,"ghost nodes involved in extrapolation. bad. \n");
        if (philevel[otherNodes[0]]<0 && philevel[otherNodes[1]]<0 && philevel[otherNodes[2]]<0)
          fprintf(stderr,"ghost nodes inevitablly involved in extrapolation. too bad. \n");

      //  break;
      }
      // check how many images found   
      if (howmany!=1) fprintf(stderr, "found %d images.\n", howmany);   
      
      // update state vector
      for (int i=0; i<dim; i++) U[iNode][i] = tempU[i];
      fprintf(stderr,"updated. U[%d] = [%f, %f, %f, %f, %f].\n", iNode, U[iNode][0], U[iNode][1],
              U[iNode][2], U[iNode][3], U[iNode][4]);
      






      delete[] nbElems;
    } 

/*  fprintf(stderr, "----------------------------------\n");
  for (int iNode=0; iNode<subDomain->numNodes(); iNode++)
    if ((!realNodeTag[iNode])&&realNodeTag0[iNode]) {
      fprintf(stderr,"found new ghost node: #%d (%d->%d) [%f, %f, %f].\n", iNode, realNodeTag0[iNode], realNodeTag[iNode],
              X[iNode][0], X[iNode][1], X[iNode][2]);
    }
*/
}
























