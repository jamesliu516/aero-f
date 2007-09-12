#include <Tet.h>

#include <Edge.h>
#include <Face.h>
#include <Vector3D.h>
#include <Vector.h>
#include <BinFileHandler.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

const double Tet::third = 1.0/3.0;
const double Tet::fourth = 1.0/4.0;
const double Tet::sixth = 1.0/6.0;
const int Tet::edgeEnd[6][2] = { {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} };
const int Tet::edgeFace[6][2] = { {0,1}, {2,0}, {1,2}, {0,3}, {3,1}, {2,3} };
const int Tet::faceDef[4][3] = { {0,1,2}, {0,3,1}, {0,2,3}, {1,3,2} };

//------------------------------------------------------------------------------

TetSet::TetSet(int value)  
{

  numTets = value;
  tets = new Tet[value];

}

//------------------------------------------------------------------------------

TetSet::~TetSet()  
{
  
  if (tets) delete[] tets;

}

//------------------------------------------------------------------------------

int TetSet::read(BinFileHandler &file, int numRanges, int (*ranges)[2], int *map)  
{

  // read in number of tets in cluster (not used)
  int numClusTets;
  file.read(&numClusTets, 1);

  // read in the offset for the first tet
  BinFileHandler::OffType start, tocStart;
  file.read(&start, 1);

  // set the location of first tet in the TOC
  tocStart = file.tell();

  int count = 0;

  // read in ranges
  for (int iRange = 0; iRange < numRanges; ++iRange) {

    // compute number of tets in range
    int nTets = ranges[iRange][1] - ranges[iRange][0] + 1;

    // seek to correct position in file thanks to the table of contents
    file.seek(tocStart + sizeof(int) * ranges[iRange][0]);
    int toc;
    file.read(&toc, 1);
    file.seek(start + sizeof(int) * toc);

    for (int i = 0; i < nTets; i++)  {
      // read in the element type (not used)
      int type, volume_id;
      file.read(&type, 1);
      // read in global tet number
      file.read(map + count, 1);
      // read in volume id (for porous)
      file.read(&volume_id, 1);
      tets[count].setVolumeID(volume_id);
      // read in tet nodes
      file.read(&tets[count][0], 4);
      count++;
    }
  }

  if (count != numTets) {
    fprintf(stderr, "*** Error: wrong number of tets read (%d instead of %d)\n",
	    count, numTets);
    exit(1);
  }

  return numClusTets;

}

//------------------------------------------------------------------------------

double Tet::computeLongestEdge(SVec<double,3> &X)
{

  double result = 0.0;
  double temp;
  Vec3D x[4] = {X[ nodeNum[0] ], X[ nodeNum[1] ], X[ nodeNum[2] ], X[ nodeNum[3] ]};

  Vec3D e[6];
  e[0] = x[0]-x[1];
  e[1] = x[0]-x[2];
  e[2] = x[0]-x[3];
  e[3] = x[1]-x[2];
  e[4] = x[1]-x[3];
  e[5] = x[2]-x[3];

  for (int i=0; i<6; i++){
      temp = sqrt(e[i][0]*e[i][0]+e[i][1]*e[i][1]+e[i][2]*e[i][2]);
      if (temp > result) result = temp;
  }
  return result;

}

//------------------------------------------------------------------------------

double Tet::computeVolume(SVec<double,3> &X)
{

  Vec3D x[4] = {X[ nodeNum[0] ], X[ nodeNum[1] ], X[ nodeNum[2] ], X[ nodeNum[3] ]};

  Vec3D v1 = x[1] - x[0];
  Vec3D v2 = x[2] - x[0];
  Vec3D v3 = x[3] - x[0];

  double volume = sixth * (v3 * (v1 ^ v2));

  return volume;

}

//------------------------------------------------------------------------------

double Tet::computeGeometricVolume(Vec3D P0, Vec3D P1, Vec3D P2, Vec3D P3)
{

  Vec3D v1 = P1-P0;
  Vec3D v2 = P2-P0;
  Vec3D v3 = P3-P0;

  double volume = sixth * (v3 * (v1 ^ v2));

  return fabs(volume);

}
//------------------------------------------------------------------------------

double Tet::computeControlVolumes(SVec<double,3> &X, Vec<double> &ctrlVol)
{

  double volume = computeVolume(X);

  for (int j=0; j<4; ++j) 
    ctrlVol[ nodeNum[j] ] += 0.25 * volume;

  return volume;

}

//------------------------------------------------------------------------------

void Tet::printInvalidElement(int numInvElem, double lscale, int i, int *nodeMap, 
			      int *tetMap, SVec<double,3> &x0, SVec<double,3> &x)
{

  int glno = tetMap[i]+1;

//  fprintf(stderr, "*** Error: negative volume for tetrahedron %d\n", glno);

  if (numInvElem > 10) return;

  char name[MAXLINE];
  sprintf(name, "tet%d.xpost", glno);

  FILE *fp = fopen(name, "w");

  if (!fp)
    fprintf(stderr, "*** Error: could not open \'%s\'\n", name);
  else {
    int j;
    fprintf(fp, "Nodes Nodes%d\n", glno);
    for (j=0; j<4; ++j) 
      fprintf(fp, "%d %e %e %e\n", j+1, lscale*x0[ nodeNum[j] ][0], 
	      lscale*x0[ nodeNum[j] ][1], lscale*x0[ nodeNum[j] ][2]);
    for (j=0; j<4; ++j) 
      fprintf(fp, "%d %e %e %e\n", j+5, lscale*x[ nodeNum[j] ][0], 
	      lscale*x[ nodeNum[j] ][1], lscale*x[ nodeNum[j] ][2]);
    fprintf(fp, "Elements Element%d using Nodes%d\n", glno, glno);
    fprintf(fp, "1 5 1 2 3 4\n");
    fprintf(fp, "Elements BadElement%d using Nodes%d\n", glno, glno);
    fprintf(fp, "1 5 5 6 7 8\n");
    fflush(fp);
    fclose(fp);
  }

  fflush(stderr);

}

//------------------------------------------------------------------------------

void Tet::numberEdges(EdgeSet &edges)
{

  edgeNum[0] = edges.find(nodeNum[ edgeEnd[0][0] ], nodeNum[ edgeEnd[0][1] ]);
  edgeNum[1] = edges.find(nodeNum[ edgeEnd[1][0] ], nodeNum[ edgeEnd[1][1] ]);
  edgeNum[2] = edges.find(nodeNum[ edgeEnd[2][0] ], nodeNum[ edgeEnd[2][1] ]);
  edgeNum[3] = edges.find(nodeNum[ edgeEnd[3][0] ], nodeNum[ edgeEnd[3][1] ]);
  edgeNum[4] = edges.find(nodeNum[ edgeEnd[4][0] ], nodeNum[ edgeEnd[4][1] ]);
  edgeNum[5] = edges.find(nodeNum[ edgeEnd[5][0] ], nodeNum[ edgeEnd[5][1] ]);

}

//------------------------------------------------------------------------------

void Tet::renumberEdges(Vec<int> &newNum)
{

  for (int l=0; l<6; ++l) 
    edgeNum[l] = newNum[ edgeNum[l] ];

}

//------------------------------------------------------------------------------

int Tet::countNodesOnBoundaries(Vec<bool> &tagNodes)
{

  int num = 0;

  if (tagNodes[ nodeNum[0] ]) ++num;
  if (tagNodes[ nodeNum[1] ]) ++num;
  if (tagNodes[ nodeNum[2] ]) ++num;
  if (tagNodes[ nodeNum[3] ]) ++num;

  return num;

}

//------------------------------------------------------------------------------

int Tet::setFaceToElementConnectivity(int i, Vec<bool> &tagNodes, 
				       MapFaces &mf, FaceSet &faces)
{

  int nswap = 0;

  for (int k=0; k<4; ++k) {

    int nn[3] = {nodeNum[ faceDef[k][0] ], 
		 nodeNum[ faceDef[k][1] ], 
		 nodeNum[ faceDef[k][2] ]};

    if (tagNodes[ nn[0] ] && tagNodes[ nn[1] ] && tagNodes[ nn[2] ]) {
      Face f;
      f.setup(2, nn);
      f.reorder();
      MapFaces::iterator it = mf.find(f);
      if (it != mf.end()) {
	int fn = (*it).second;
	if ((faces[fn][0] == nn[0] && faces[fn][1] == nn[1] && faces[fn][2] == nn[2]) ||
	    (faces[fn][0] == nn[1] && faces[fn][1] == nn[2] && faces[fn][2] == nn[0]) ||
	    (faces[fn][0] == nn[2] && faces[fn][1] == nn[0] && faces[fn][2] == nn[1]))
	  faces[fn].setElementNumber(i, 0);
	else {
	  ++nswap;
	  faces[fn].setElementNumber(i, nn);
	}
      }
    }
  }

  return nswap;

}

//------------------------------------------------------------------------------
// define the calculation of edge normal for each tetrahedron
/*
              
       g__________ 0
       /\        /
      /  \      /
     /    \    /
    /      \  /
  1/________\/
             e

   g = CG of tetrahedron
   e = mid-edge
   0 = CG of face 0
   1 = CG of face 1
	     
*/
/*
@ARTICLE{koobus-farhat-99a,
  author = "Koobus, B. and Farhat, C.",
  title = "Second-order time-accurate and geometrically conservative implicit
  schemes for flow computations on unstructured dynamic meshes",
  journal = cmame,
  year = 1999,
  volume = 170,
  pages = "103--130",
}
@ARTICLE{geuzaine-grandmont-farhat-03,
  author = "Geuzaine, P. and Grandmont, C. and Farhat, C.",
  title = "Design and analysis of {ALE} schemes with provable second-order
           time-accuracy for inviscid and viscous flow simulations",
  journal = jcp,
  volume = 191,
  number = 1,
  pages = "206--227",
  year = 2003,
} 
*/
void Tet::computeEdgeNormalsGCL1(SVec<double,3> &Xn, SVec<double,3> &Xnp1, 
				 SVec<double,3> &Xdot, Vec<Vec3D> &edgeNorm, 
				 Vec<double> &edgeNormVel)
{

  static const double c0 = 13.0/36.0, c1 = 5.0/36.0;
  static const int edgeOpEnd[6][2] = { {2,3}, {3,1}, {1,2}, {0,3}, {2,0}, {0,1} };

  Vec3D x_n[4] = {Xn[ nodeNum[0] ], Xn[ nodeNum[1] ], 
		  Xn[ nodeNum[2] ], Xn[ nodeNum[3] ]};
  Vec3D x_np1[4] = {Xnp1[ nodeNum[0] ], Xnp1[ nodeNum[1] ], 
		    Xnp1[ nodeNum[2] ], Xnp1[ nodeNum[3] ]};
  Vec3D xdot[4] = {Xdot[ nodeNum[0] ], Xdot[ nodeNum[1] ], 
		   Xdot[ nodeNum[2] ], Xdot[ nodeNum[3] ]};

  Vec3D f_n[4];
  f_n[0] = third * (x_n[0] + x_n[1] + x_n[2]);
  f_n[1] = third * (x_n[0] + x_n[1] + x_n[3]);
  f_n[2] = third * (x_n[0] + x_n[2] + x_n[3]);
  f_n[3] = third * (x_n[1] + x_n[2] + x_n[3]);

  Vec3D g_n = 0.25 * (x_n[0] + x_n[1] + x_n[2] + x_n[3]);

  Vec3D f_np1[4];
  f_np1[0] = third * (x_np1[0] + x_np1[1] + x_np1[2]);
  f_np1[1] = third * (x_np1[0] + x_np1[1] + x_np1[3]);
  f_np1[2] = third * (x_np1[0] + x_np1[2] + x_np1[3]);
  f_np1[3] = third * (x_np1[1] + x_np1[2] + x_np1[3]);

  Vec3D g_np1 = 0.25 * (x_np1[0] + x_np1[1] + x_np1[2] + x_np1[3]);

  for (int l=0; l<6; ++l) {
    Vec3D xg0_n = f_n[ edgeFace[l][0] ] - g_n;
    Vec3D xg1_n = f_n[ edgeFace[l][1] ] - g_n;
    Vec3D xg0_np1 = f_np1[ edgeFace[l][0] ] - g_np1;
    Vec3D xg1_np1 = f_np1[ edgeFace[l][1] ] - g_np1;

    Vec3D n = 0.5 * ((xg1_np1 ^ xg0_np1) + (xg1_n ^ xg0_n)) + 
      0.25 * ((xg1_np1 ^ xg0_n) + (xg1_n ^ xg0_np1));

    double ndot = ( c0 * (xdot[ edgeEnd[l][0] ] + xdot[ edgeEnd[l][1] ]) + 
		    c1 * (xdot[ edgeOpEnd[l][0] ] + xdot[ edgeOpEnd[l][1] ]) ) * n;

    if (nodeNum[ edgeEnd[l][0] ] < nodeNum[ edgeEnd[l][1] ]) {
      edgeNorm[ edgeNum[l] ] += n;
      edgeNormVel[ edgeNum[l] ] += ndot;
    }
    else {
      edgeNorm[ edgeNum[l] ] -= n;
      edgeNormVel[ edgeNum[l] ] -= ndot;
    }
  }

}

//------------------------------------------------------------------------------

void Tet::computeEdgeNormalsEZGCL1(double oodt, SVec<double,3> &Xn, SVec<double,3> &Xnp1, 
				   Vec<Vec3D> &edgeNorm, Vec<double> &edgeNormVel)
{

  Vec3D x_n[4] = {Xn[ nodeNum[0] ], Xn[ nodeNum[1] ], 
		  Xn[ nodeNum[2] ], Xn[ nodeNum[3] ]};
  Vec3D x_np1[4] = {Xnp1[ nodeNum[0] ], Xnp1[ nodeNum[1] ], 
		    Xnp1[ nodeNum[2] ], Xnp1[ nodeNum[3] ]};

  Vec3D e_n[6], f_n[4], g_n;

  e_n[0] = 0.5 * (x_n[0] + x_n[1]);
  e_n[1] = 0.5 * (x_n[0] + x_n[2]);
  e_n[2] = 0.5 * (x_n[0] + x_n[3]);
  e_n[3] = 0.5 * (x_n[1] + x_n[2]);
  e_n[4] = 0.5 * (x_n[1] + x_n[3]);
  e_n[5] = 0.5 * (x_n[2] + x_n[3]);

  f_n[0] = third * (x_n[0] + x_n[1] + x_n[2]);
  f_n[1] = third * (x_n[0] + x_n[1] + x_n[3]);
  f_n[2] = third * (x_n[0] + x_n[2] + x_n[3]);
  f_n[3] = third * (x_n[1] + x_n[2] + x_n[3]);

  g_n = 0.5 * (e_n[0] + e_n[5]);

  Vec3D e_np1[6], f_np1[4], g_np1;

  e_np1[0] = 0.5 * (x_np1[0] + x_np1[1]);
  e_np1[1] = 0.5 * (x_np1[0] + x_np1[2]);
  e_np1[2] = 0.5 * (x_np1[0] + x_np1[3]);
  e_np1[3] = 0.5 * (x_np1[1] + x_np1[2]);
  e_np1[4] = 0.5 * (x_np1[1] + x_np1[3]);
  e_np1[5] = 0.5 * (x_np1[2] + x_np1[3]);

  f_np1[0] = third * (x_np1[0] + x_np1[1] + x_np1[2]);
  f_np1[1] = third * (x_np1[0] + x_np1[1] + x_np1[3]);
  f_np1[2] = third * (x_np1[0] + x_np1[2] + x_np1[3]);
  f_np1[3] = third * (x_np1[1] + x_np1[2] + x_np1[3]);

  g_np1 = 0.5 * (e_np1[0] + e_np1[5]);

  for (int l=0; l<6; ++l) {
    Vec3D e0_np1 = f_np1[ edgeFace[l][0] ] - e_np1[l];
    Vec3D e1_np1 = f_np1[ edgeFace[l][1] ] - e_np1[l];
    Vec3D eg_np1 = g_np1 - e_np1[l];
    /*EZ1
    Vec3D e0_np1 = f_n[ edgeFace[l][0] ] - e_n[l];
    Vec3D e1_np1 = f_n[ edgeFace[l][1] ] - e_n[l];
    Vec3D eg_np1 = g_n - e_n[l];
    */
    Vec3D n0 = 0.5 * (e0_np1 ^ eg_np1);
    Vec3D n1 = 0.5 * (eg_np1 ^ e1_np1);
    Vec3D n = n0 + n1;

    double vol0 = Face::computeVolume(e_n[l], f_n[ edgeFace[l][0] ], g_n,
				      e_np1[l], f_np1[ edgeFace[l][0] ], g_np1);
    double vol1 = Face::computeVolume(e_n[l], g_n, f_n[ edgeFace[l][1] ],
				      e_np1[l], g_np1, f_np1[ edgeFace[l][1] ]);
    double ndot = oodt * (vol0 + vol1);

    if (nodeNum[ edgeEnd[l][0] ] < nodeNum[ edgeEnd[l][1] ]) {
      edgeNorm[ edgeNum[l] ] += n;
      edgeNormVel[ edgeNum[l] ] += ndot;
    }
    else {
      edgeNorm[ edgeNum[l] ] -= n;
      edgeNormVel[ edgeNum[l] ] -= ndot;
    }
  }

}

//------------------------------------------------------------------------------
/*
void Tet::computeEdgeNormalsLZGCL1(SVec<double,3> &Xn, SVec<double,3> &Xnp1, 
				   SVec<double,3> &Xdot, Vec<Vec3D> &edgeNorm, 
				   Vec<double> &edgeNormVel)
{

  static double twelfth = 1.0/12.0;

  Vec3D x_n[4] = {Xn[ nodeNum[0] ], Xn[ nodeNum[1] ], 
		  Xn[ nodeNum[2] ], Xn[ nodeNum[3] ]};
  Vec3D x_np1[4] = {Xnp1[ nodeNum[0] ], Xnp1[ nodeNum[1] ], 
		    Xnp1[ nodeNum[2] ], Xnp1[ nodeNum[3] ]};
  Vec3D xdot[4] = {Xdot[ nodeNum[0] ], Xdot[ nodeNum[1] ], 
		   Xdot[ nodeNum[2] ], Xdot[ nodeNum[3] ]};

  Vec3D midEdge_n[6], cgFace_n[4], cgTet_n;

  midEdge_n[0] = 0.5 * (x_n[0] + x_n[1]);
  midEdge_n[1] = 0.5 * (x_n[0] + x_n[2]);
  midEdge_n[2] = 0.5 * (x_n[0] + x_n[3]);
  midEdge_n[3] = 0.5 * (x_n[1] + x_n[2]);
  midEdge_n[4] = 0.5 * (x_n[1] + x_n[3]);
  midEdge_n[5] = 0.5 * (x_n[2] + x_n[3]);

  cgFace_n[0] = third * (x_n[0] + x_n[1] + x_n[2]);
  cgFace_n[1] = third * (x_n[0] + x_n[1] + x_n[3]);
  cgFace_n[2] = third * (x_n[0] + x_n[2] + x_n[3]);
  cgFace_n[3] = third * (x_n[1] + x_n[2] + x_n[3]);

  cgTet_n = 0.5 * (midEdge_n[0] + midEdge_n[5]);

  Vec3D midEdge_np1[6], cgFace_np1[4], cgTet_np1;

  midEdge_np1[0] = 0.5 * (x_np1[0] + x_np1[1]);
  midEdge_np1[1] = 0.5 * (x_np1[0] + x_np1[2]);
  midEdge_np1[2] = 0.5 * (x_np1[0] + x_np1[3]);
  midEdge_np1[3] = 0.5 * (x_np1[1] + x_np1[2]);
  midEdge_np1[4] = 0.5 * (x_np1[1] + x_np1[3]);
  midEdge_np1[5] = 0.5 * (x_np1[2] + x_np1[3]);

  cgFace_np1[0] = third * (x_np1[0] + x_np1[1] + x_np1[2]);
  cgFace_np1[1] = third * (x_np1[0] + x_np1[1] + x_np1[3]);
  cgFace_np1[2] = third * (x_np1[0] + x_np1[2] + x_np1[3]);
  cgFace_np1[3] = third * (x_np1[1] + x_np1[2] + x_np1[3]);

  cgTet_np1 = 0.5 * (midEdge_np1[0] + midEdge_np1[5]);

  Vec3D midEdgeVel[6], cgFaceVel[4], cgTetVel;

  midEdgeVel[0] = 0.5 * (xdot[0] + xdot[1]);
  midEdgeVel[1] = 0.5 * (xdot[0] + xdot[2]);
  midEdgeVel[2] = 0.5 * (xdot[0] + xdot[3]);
  midEdgeVel[3] = 0.5 * (xdot[1] + xdot[2]);
  midEdgeVel[4] = 0.5 * (xdot[1] + xdot[3]);
  midEdgeVel[5] = 0.5 * (xdot[2] + xdot[3]);

  cgFaceVel[0] = third * (xdot[0] + xdot[1] + xdot[2]);
  cgFaceVel[1] = third * (xdot[0] + xdot[1] + xdot[3]);
  cgFaceVel[2] = third * (xdot[0] + xdot[2] + xdot[3]);
  cgFaceVel[3] = third * (xdot[1] + xdot[2] + xdot[3]);

  cgTetVel = 0.5 * (midEdgeVel[0] + midEdgeVel[5]);

  for (int l=0; l<6; ++l) {

    Vec3D e0_n = cgFace_n[ edgeFace[l][0] ] - midEdge_n[l];
    Vec3D e1_n = cgFace_n[ edgeFace[l][1] ] - midEdge_n[l];
    Vec3D eg_n = cgTet_n - midEdge_n[l];

    Vec3D e0_np1 = cgFace_np1[ edgeFace[l][0] ] - midEdge_np1[l];
    Vec3D e1_np1 = cgFace_np1[ edgeFace[l][1] ] - midEdge_np1[l];
    Vec3D eg_np1 = cgTet_np1 - midEdge_np1[l];

    Vec3D vel0 = third * (midEdgeVel[l] + cgTetVel + cgFaceVel[ edgeFace[l][0] ]);
    Vec3D vel1 = third * (midEdgeVel[l] + cgTetVel + cgFaceVel[ edgeFace[l][1] ]);
    
    Vec3D n0 = twelfth * (((2.0*e0_np1 + e0_n) ^ eg_np1) + ((2.0*e0_n + e0_np1) ^ eg_n));
    Vec3D n1 = twelfth * (((2.0*eg_np1 + eg_n) ^ e1_np1) + ((2.0*eg_n + eg_np1) ^ e1_n));

    Vec3D n = n0 + n1;
    double ndot = vel0 * n0 + vel1 * n1;

    if (nodeNum[ edgeEnd[l][0] ] < nodeNum[ edgeEnd[l][1] ]) {
      edgeNorm[ edgeNum[l] ] += n;
      edgeNormVel[ edgeNum[l] ] += ndot;
    }
    else {
      edgeNorm[ edgeNum[l] ] -= n;
      edgeNormVel[ edgeNum[l] ] -= ndot;
    }

  }

}
*/
//------------------------------------------------------------------------------

void Tet::computeWeightsGalerkin(SVec<double,3> &X, SVec<double,3> &wii,
				 SVec<double,3> &wij, SVec<double,3> &wji)
{

  int i, j, k;

  double dp1dxj[4][3];

  double vol = computeGradientP1Function(X, dp1dxj);

  for (k=0; k<4; ++k) {

    dp1dxj[k][0] *= vol;
    dp1dxj[k][1] *= vol;
    dp1dxj[k][2] *= vol;

    wii[ nodeNum[k] ][0] += dp1dxj[k][0];
    wii[ nodeNum[k] ][1] += dp1dxj[k][1];
    wii[ nodeNum[k] ][2] += dp1dxj[k][2];

  }

  for (k=0; k<6; ++k) {

    if (nodeNum[ edgeEnd[k][0] ] < nodeNum[ edgeEnd[k][1] ]) {
      i = edgeEnd[k][0];
      j = edgeEnd[k][1];
    } 
    else {
      i = edgeEnd[k][1];
      j = edgeEnd[k][0];
    }

    wji[ edgeNum[k] ][0] += dp1dxj[i][0];
    wij[ edgeNum[k] ][0] += dp1dxj[j][0];

    wji[ edgeNum[k] ][1] += dp1dxj[i][1];
    wij[ edgeNum[k] ][1] += dp1dxj[j][1];

    wji[ edgeNum[k] ][2] += dp1dxj[i][2];
    wij[ edgeNum[k] ][2] += dp1dxj[j][2];

  }

}

//------------------------------------------------------------------------------

void Tet::computeEdgeWeightsGalerkin(SVec<double,3> &X, SVec<double,9> &M)
{

  static int faces[4][3] = { {0,1,3}, {0,2,1}, {1,2,3}, {0,3,2} };
  static int edges[6][2] = { {3,2}, {0,2}, {1,2}, {0,3}, {1,3}, {1,0} };

  Vec3D x[4] = {X[ nodeNum[0] ], X[ nodeNum[1] ], X[ nodeNum[2] ], X[ nodeNum[3] ]};

  double invvol = 1.0 / computeVolume(X);

  for (int l=0; l<6; ++l) {

    Vec3D n0 = 0.5 * ((x[ faces[ edges[l][0] ][1] ] - x[ faces[ edges[l][0] ][0] ]) ^ 
		      (x[ faces[ edges[l][0] ][2] ] - x[ faces[ edges[l][0] ][0] ]));
    Vec3D n1 = 0.5 * ((x[ faces[ edges[l][1] ][1] ] - x[ faces[ edges[l][1] ][0] ]) ^ 
		      (x[ faces[ edges[l][1] ][2] ] - x[ faces[ edges[l][1] ][0] ]));

    if (nodeNum[ edgeEnd[l][0] ] > nodeNum[ edgeEnd[l][1] ]) {
      Vec3D tmp = n1;
      n1 = n0;
      n0 = tmp;
    }

    M[ edgeNum[l] ][0] += invvol * n1[0] * n0[0];
    M[ edgeNum[l] ][1] += invvol * n1[0] * n0[1];
    M[ edgeNum[l] ][2] += invvol * n1[0] * n0[2];
    M[ edgeNum[l] ][3] += invvol * n1[1] * n0[1];
    M[ edgeNum[l] ][4] += invvol * n1[1] * n0[2];
    M[ edgeNum[l] ][5] += invvol * n1[2] * n0[2];
    
    M[ edgeNum[l] ][6] += invvol * n1[1] * n0[0];
    M[ edgeNum[l] ][7] += invvol * n1[2] * n0[0];
    M[ edgeNum[l] ][8] += invvol * n1[2] * n0[1];

  }

}

//------------------------------------------------------------------------------

void Tet::computeStiffAndForce(double force[4][3], double K[12][12],
			       SVec<double,3> &X, SVec<double,3> &nodes, double volStiff)
{

  // X is the current position of nodes
  // nodes is the reference position of nodes

  int i, j, k;
  double nGrad[4][3];

  // cast to simplify loops:
  double *f = reinterpret_cast<double *>(force);

  // compute dN_i/dX_j, also obtain dOmega 
  // (actually 1/4th of it since we have a factor 2 on e and s

  double realVol = computeGradientP1Function(nodes, nGrad);
  
  // Scaling of this stiffness for aeroelastic reasons:
  //dOmega = pow(dOmega, 2.0/3.0);
  // Remove volume scaling => This gives small elements more stiffness
  double dOmega = 1.0;

  double V[12];
  double VV[12][12];

  // now get F_ij = dPhi_i/dX_j = x^k_i dN_k/dX_j
  double F[3][3];
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      F[i][j] = X[nodeNum[0]][i]*nGrad[0][j] +
  	        X[nodeNum[1]][i]*nGrad[1][j] +
                X[nodeNum[2]][i]*nGrad[2][j] +
                X[nodeNum[3]][i]*nGrad[3][j];

  // compute e_ij = F_ki Fkj - delta_ij
  // This is really 2*e, but that means we simply have a factor 4 in
  //all our results
  double e_11 = F[0][0]*F[0][0]+F[1][0]*F[1][0]+F[2][0]*F[2][0] - 1.0;
  double e_22 = F[0][1]*F[0][1]+F[1][1]*F[1][1]+F[2][1]*F[2][1] - 1.0;
  double e_33 = F[0][2]*F[0][2]+F[1][2]*F[1][2]+F[2][2]*F[2][2] - 1.0;
  double e_12 = F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1];
  double e_13 = F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2];
  double e_23 = F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2];
  double sigma[6];
  double nu = 0.33;
  double E = 1.0;
  // E2 == E/(1+nu)   E1 = lambda+2*G ==  E*(1-nu)/((1+nu)(1-2nu))

  double E2 = E*nu/((1+nu)*(1-2*nu));
  double G2 = E/(1+nu);
  double E1 = E2+E/(1+nu);
  sigma[0] = E1*e_11+E2*(e_22+e_33);
  sigma[1] = E1*e_22+E2*(e_11+e_33);
  sigma[2] = E1*e_33+E2*(e_11+e_22);
  sigma[3] = 2.0*G2*e_12;
  sigma[4] = 2.0*G2*e_13;
  sigma[5] = 2.0*G2*e_23;

  // Compute de_ij/dUl for the symmetric part
  // First we get df_ij/dUl in a very compact form.
  // df_ij/dUl = dN_p/dX_j delta_iq; with p = int(l/3)+1 and l-1=q-1 mod(3)
  // this means that df_ij/dUl is already contained in dN_k/dX_j
  double dedU[12][6];
  for (i = 0; i < 4; ++i)
    for (j = 0; j < 3; ++j) {
      dedU[3*i+j][0] = 2.0*nGrad[i][0]*F[j][0];
      dedU[3*i+j][1] = 2.0*nGrad[i][1]*F[j][1];
      dedU[3*i+j][2] = 2.0*nGrad[i][2]*F[j][2];
      dedU[3*i+j][3] = nGrad[i][0]*F[j][1] + nGrad[i][1]*F[j][0];
      dedU[3*i+j][4] = nGrad[i][0]*F[j][2] + nGrad[i][2]*F[j][0];
      dedU[3*i+j][5] = nGrad[i][1]*F[j][2] + nGrad[i][2]*F[j][1];
    }

  // Get the force:

  for (i = 0; i < 12; ++i)  {
    f[i] = dOmega*( dedU[i][0]*sigma[0] + dedU[i][1]*sigma[1] +
	            dedU[i][2]*sigma[2] + dedU[i][3]*sigma[3] +
		    dedU[i][4]*sigma[4] + dedU[i][5]*sigma[5]);
  }

  // now get ds_ij/dUl
  double dsdU[12][6];
  for (i = 0; i < 12; ++i) {
    dsdU[i][0] = E1*dedU[i][0]+E2*(dedU[i][1]+dedU[i][2]);
    dsdU[i][1] = E1*dedU[i][1]+E2*(dedU[i][0]+dedU[i][2]);
    dsdU[i][2] = E1*dedU[i][2]+E2*(dedU[i][0]+dedU[i][1]);
    // the shear terms are doubled to have the full effect
    dsdU[i][3] = 2.0*G2*dedU[i][3];
    dsdU[i][4] = 2.0*G2*dedU[i][4];
    dsdU[i][5] = 2.0*G2*dedU[i][5];
  }

  // multiply modified dsdU by dedU Only do the symmetric part
  for (i = 0; i < 12; ++i)
    for (j = 0; j <= i; ++j)
      K[j][i]  = dsdU[i][0]*dedU[j][0] + dsdU[i][1]*dedU[j][1] +
                 dsdU[i][2]*dedU[j][2] + dsdU[i][3]*dedU[j][3] +
                 dsdU[i][4]*dedU[j][4] + dsdU[i][5]*dedU[j][5];

  // add s*d2e/dU_idUj (symmetric part only)
  for (i = 0; i < 4; ++i)
    for (k = 0; k <= i; ++k)
      for (j = 0; j < 3; ++j)
        K[3*k+j][3*i+j] +=
                sigma[0]*(2*nGrad[i][0]*nGrad[k][0]) +
                sigma[1]*(2*nGrad[i][1]*nGrad[k][1]) +
                sigma[2]*(2*nGrad[i][2]*nGrad[k][2]) +
                sigma[3]*(nGrad[i][0]*nGrad[k][1]+nGrad[i][1]*nGrad[k][0]) +
                sigma[4]*(nGrad[i][0]*nGrad[k][2]+nGrad[i][2]*nGrad[k][0]) +
                sigma[5]*(nGrad[i][1]*nGrad[k][2]+nGrad[i][2]*nGrad[k][1]);
		
  if (volStiff > 0.0)  {
    
    // Compute Volume
    static double sixth = 1.0/6.0;

    // compute factor of sixth divided by reference volume
    double invOmega = 1.0 / realVol;
 
    Vec3D x[4] = {X[ nodeNum[0] ], X[ nodeNum[1] ], X[ nodeNum[2] ], X[ nodeNum[3] ]};

    Vec3D v1 = x[3] - x[1];
    Vec3D v2 = x[2] - x[1];
    Vec3D v3 = x[0] - x[1];

    double volume = sixth * (v3 * (v1 ^ v2));
    if (volume < 0.0) {
      fprintf(stderr, "*** Error: negative jacobian\n");
      exit(-1);
    }

    // compute volume derivatives 
    // dV/dX0 OK
    Vec3D cross = (v1 ^ v2);
    V[0] = sixth * cross[0];
    V[1] = sixth * cross[1];
    V[2] = sixth * cross[2];
  
    // dV/dX1
    v1 = x[2] - x[0];
    v2 = x[3] - x[0];
    cross = (v1 ^ v2);
    V[3] = sixth * (cross[0]);
    V[4] = sixth * (cross[1]);
    V[5] = sixth * (cross[2]);
  
    // dV/dX2
    v1 = x[3] - x[0];
    v2 = x[1] - x[0];
    cross = (v1 ^ v2);
    V[6] = sixth * (cross[0]);
    V[7] = sixth * (cross[1]);
    V[8] = sixth * (cross[2]);
 
    // dV/dX3
    v1 = x[1] - x[0];
    v2 = x[2] - x[0];
    cross = (v1 ^ v2);
    V[9] = sixth * (cross[0]);
    V[10] = sixth * (cross[1]);
    V[11] = sixth * (cross[2]);
    
    // compute force contribution from volume energy function
    // dW(V)/dX == (2*(V/V0 - 1)*(1/V0)*(dV/dX)) / ((V/V0)*(V/V0)*(V/V0))
    double volRatio = volume * invOmega;

    double coef = volStiff*(dOmega*2.0 *(volRatio - 1.0) * invOmega) / (volRatio*volRatio*volRatio);

    for (i = 0; i < 12; i++)
      f[i] += coef * V[i];

    // compute 2nd derivatives of V
    // just to make sure
    for (i = 0; i< 12;++i)
      for (j=0;j<12;j++)
        VV[i][j] = 0.0;

    v1 = sixth * (x[2]-x[3]); 
    VV[0][3] = 0;
    VV[0][4] = -v1[2];
    VV[0][5] = v1[1];
    VV[1][3] = v1[2];
    VV[1][4] = 0;
    VV[1][5] = -v1[0];
    VV[2][3] = -v1[1];
    VV[2][4] = v1[0];
    VV[2][5] = 0;

    v1 = sixth * (x[3] - x[1]);
    VV[0][6] = 0;
    VV[0][7] = -v1[2];
    VV[0][8] = v1[1];
    VV[1][6] = v1[2];
    VV[1][7] = 0;
    VV[1][8] = -v1[0];
    VV[2][6] = -v1[1];
    VV[2][7] = v1[0];
    VV[2][8] = 0;

    v1 = sixth * (x[1] - x[2]);
    VV[0][9] = 0;
    VV[0][10] =-v1[2];
    VV[0][11] =v1[1];
    VV[1][9] = v1[2];
    VV[1][10] =0;
    VV[1][11] =-v1[0];
    VV[2][9] = -v1[1];
    VV[2][10] =v1[0];
    VV[2][11] =0;

    v1 = sixth * (x[0] - x[3]);
    VV[3][6] = 0;
    VV[3][7] = -v1[2];
    VV[3][8] = v1[1];
    VV[4][6] = v1[2];
    VV[4][7] = 0;
    VV[4][8] = -v1[0];
    VV[5][6] = -v1[1];
    VV[5][7] = v1[0];
    VV[5][8] = 0;

    v1 = sixth * (x[2] - x[0]);
    VV[3][9] = 0;
    VV[3][10] =-v1[2];
    VV[3][11] =v1[1];
    VV[4][9] = v1[2];
    VV[4][10] =0;
    VV[4][11] =-v1[0];
    VV[5][9] = -v1[1];
    VV[5][10] =v1[0];
    VV[5][11] =0;

    v1 = sixth * (x[0] - x[1]);
    VV[6][9] = 0;
    VV[6][10] =-v1[2];
    VV[6][11] =v1[1];
    VV[7][9] = v1[2];
    VV[7][10] =0;
    VV[7][11] =-v1[0];
    VV[8][9] = -v1[1];
    VV[8][10] =v1[0];
    VV[8][11] =0;

    // Compute symmetric terms of K_ij 

    for (i = 0; i < 4; ++i)
      for (k = 0; k <= i; ++k)
        for (j = 0; j < 3; ++j)
    	  for (int l = 0; l < 3; ++l)
            K[3*k+j][3*i+l] += volStiff*((2.0 * invOmega) / (volRatio*volRatio*volRatio*volRatio))
                            * ( (3.0 - 2.0*volRatio)*invOmega*V[3*k+j]*V[3*i+l]
                            +   (volRatio - 1.0)*volRatio* VV[3*k+j][3*i+l] );

  }  

  // Symmetrize and multiply by the volume
  for (i = 0; i < 12; ++i)
    for (j = 0; j <= i; ++j)
      K[i][j] = (K[j][i] *= dOmega);

}

//------------------------------------------------------------------------------

void Tet::computeStiffAndForceLIN(double force[4][3], double K[12][12],
				  SVec<double,3> &X, SVec<double,3> &nodes)

{

  int i, j;

  // cast to simplify loops:
  double *f = reinterpret_cast<double *>(force);

  // Force is zero for linear element
  for(i = 0; i < 12; ++i)
    f[i] = 0.0;


  // compute dN_i/dX_j, also obtain dOmega (actually 1/4th of it since we have a factor 2 on e and s
  double nGrad[4][3];
  //double dOmega = fourth * computeGradientP1Function(nodes, nGrad);
  computeGradientP1Function(nodes, nGrad);

  // Scaling of this stiffness for aeroelastic reasons:
  //dOmega = pow(dOmega, 2.0/3.0);
  //dOmega = 1.0;

/*
  //  compute E, nu according to lambda = 2/3 mu = 1/Vol
  // get current volume
  static double sixth = 1.0/6.0;

  Vec3D x[4] = {X[ nodeNum[0] ], X[ nodeNum[1] ], X[ nodeNum[2] ], X[ nodeNum[3] ]};

  Vec3D v1 = x[1] - x[0];
  Vec3D v2 = x[2] - x[0];
  Vec3D v3 = x[3] - x[0];

  double vol = sixth * (v3 * (v1 ^ v2));
  //double nu = 0.2;
  //double E = 2.4 / vol;
*/

  // choose lame constants such that lam/mu = 1
  double nu = .25;
  double E = 2.5;

  // E2 == E/(1+nu)   E1 = lambda+2*G ==  E*(1-nu)/((1+nu)(1-2nu))
  double E2=E*nu/((1+nu)*(1-2*nu));
  double G2 = E/(1+nu);
  double E1 = E2+E/(1+nu);

  // Create Identity Matrix
  double F[3][3];
  F[0][0] = F[1][1] = F[2][2] = 1.0;
  F[0][1] = F[0][2] = F[1][0] = F[1][2] = F[2][0] = F[2][1] = 0.0;

  // Compute de_ij/dUl for the symmetric part
  // First we get df_ij/dUl in a very compact form.
  // df_ij/dUl = dN_p/dX_j delta_iq; with p = int(l/3)+1 and l-1=q-1 mod(3)
  // this means that df_ij/dUl is already contained in dN_k/dX_j
  double dedU[12][6];
  for(i = 0; i < 4; ++i)
    for(j = 0; j < 3; ++j) {
      dedU[3*i+j][0] = nGrad[i][0]*F[j][0];
      dedU[3*i+j][1] = nGrad[i][1]*F[j][1];
      dedU[3*i+j][2] = nGrad[i][2]*F[j][2];
      dedU[3*i+j][3] = .5*(nGrad[i][0]*F[j][1] + nGrad[i][1]*F[j][0]);
      dedU[3*i+j][4] = .5*(nGrad[i][0]*F[j][2] + nGrad[i][2]*F[j][0]);
      dedU[3*i+j][5] = .5*(nGrad[i][1]*F[j][2] + nGrad[i][2]*F[j][1]);
    }

  // now get ds_ij/dUl
  double dsdU[12][6];
  for(i = 0; i < 12; ++i) {
    dsdU[i][0] = E1*dedU[i][0]+E2*(dedU[i][1]+dedU[i][2]);
    dsdU[i][1] = E1*dedU[i][1]+E2*(dedU[i][0]+dedU[i][2]);
    dsdU[i][2] = E1*dedU[i][2]+E2*(dedU[i][0]+dedU[i][1]);
    // the shear terms are doubled to have the full effect
    dsdU[i][3] = 2*G2*dedU[i][3];
    dsdU[i][4] = 2*G2*dedU[i][4];
    dsdU[i][5] = 2*G2*dedU[i][5];
  }

  // multiply modified dsdU by dedU Only do the symmetric part
  for(i = 0; i < 12; ++i)
    for(j = 0; j <= i; ++j)
      K[j][i] =  dsdU[i][0]*dedU[j][0] + dsdU[i][1]*dedU[j][1]
	      +  dsdU[i][2]*dedU[j][2] + dsdU[i][3]*dedU[j][3]
              +  dsdU[i][4]*dedU[j][4] + dsdU[i][5]*dedU[j][5];

  // Symmetrize and multiply by the volume
  for(i = 0; i < 12; ++i)
    for(j = 0; j <= i; ++j)
      K[i][j] = K[j][i];
      //K[i][j] = (K[j][i] *= dOmega);
}

//------------------------------------------------------------------------------


void TetSet::computeDynamicLESTerm(DynamicLESTerm *dles, SVec<double,2> &CsDeltaSq, SVec<double,3> &X,
                                   Vec<double> &Cs, Vec<double> &VolSum)
{

 for (int i=0; i<numTets; ++i)
    tets[i].computeDynamicLESTerm(dles, CsDeltaSq, X, Cs, VolSum);

}


//------------------------------------------------------------------------------


void Tet::computeDynamicLESTerm(DynamicLESTerm *dles, SVec<double,2> &CsDeltaSq,
                                SVec<double,3> &X, Vec<double> &Cs, Vec<double> &VolSum)
{


  double dp1dxj[4][3];
  double vol = computeGradientP1Function(X, dp1dxj);

  double delta = 0.0;
  dles->computeDelta(vol, X, nodeNum, delta);


  for(int i=0; i<4; ++i) {
    if(isnan(sqrt( CsDeltaSq[nodeNum[i]][0] ))) Cs[nodeNum[i]] = CsDeltaSq[nodeNum[i]][0];
    else  Cs[nodeNum[i]] = sqrt(CsDeltaSq[nodeNum[i]][0]/VolSum[nodeNum[i]])/delta;
  }


}


//------------------------------------------------------------------------------

void Tet::computePsiResidualSubTet(double psi[4], double phi[4],
				   Vec3D A, Vec3D B, Vec3D C, Vec3D D,
				   double locdphi[4], double locw[4],
 				   double locbeta[4], bool debug)
{

  //A,B,C,D are supposed to be the nodes of an oriented tet
  //cf Barth and Sethian for method (JCP, 1998)
  // input notation is like our code
  // implementation notation is same as paper by Barth & Sethian
  for(int i=0; i<4; i++){
    locdphi[i] = 0.0;
    locw[i]    = 0.0;
    locbeta[i] = 0.0;
  }

  //compute gradient of shape functions
  double dp1dxj[4][3];
  double vol = computeGradientP1Function(A,B,C,D, dp1dxj);
  if(vol<= 0.0) fprintf(stdout, "vol = %e\n", vol);
  assert(vol>0);

  double n[4][3];
  for (int i=0; i<4; i++)
    for (int j=0; j<3; j++)
      n[i][j] = vol*dp1dxj[i][j];

  //compute of gradient of Psi in tet and its norm
  double grad[3];
  for (int i=0; i<3; i++){
    grad[i] = 0.0;
    for (int j=0; j<4; j++)
      grad[i] += dp1dxj[j][i]*psi[j];
  }
  double oonormg = sqrt(grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);
  assert(oonormg>0.0);
  oonormg = 1.0/oonormg;

  //phitet tells sign of phi in tet, 
  //it is averaged for cases when some nodes are phi=0.0
  double phitet = 0.25*(phi[0]+phi[1]+phi[2]+phi[3]);
  if(debug){
    fprintf(stdout, "gradpsi = %e %e %e\n", grad[0],grad[1],grad[2]);
    fprintf(stdout, "phitet = %e\n", phitet);
  }

  double Fadv, Fsrc;
  if(phitet>0.0)
    Fadv = 1.0;
  else
    Fadv = -1.0;
  Fsrc = Fadv;

  //compute K for each vertex of the tet
  double K[4], Kp[4], Km[4];
  for (int i=0; i<4; i++){
    K[i] = Fadv*(grad[0]*n[i][0]+grad[1]*n[i][1]+grad[2]*n[i][2])*oonormg;
    Kp[i] = max(K[i],0.0);
    Km[i] = min(K[i],0.0);
    if(debug){
      fprintf(stdout, "K = %e - Kp = %e - Km = %e\n", K[i],Kp[i],Km[i]);
    }
  }
  double oosumKm = Km[0]+Km[1]+Km[2]+Km[3];
  assert(oosumKm<0.0);
  oosumKm = 1.0/oosumKm;

  //compute deltaphi for tet
  double dphi = K[0]*psi[0]+K[1]*psi[1]+K[2]*psi[2]+K[3]*psi[3];

  //compute deltaphi for each node
  double dphii[4];
  for (int i=0; i<4; i++){
    dphii[i] = 0.0;
    for (int j=0; j<4; j++)
      dphii[i] += Km[j]*(psi[i]-psi[j]);
    dphii[i] *= (Kp[i]*oosumKm);
  }

  //compute alpha coefficients
  double alphatot = max(0.0,dphii[0]/dphi)+max(0.0,dphii[1]/dphi)+
                    max(0.0,dphii[2]/dphi)+max(0.0,dphii[3]/dphi);
  assert(alphatot!=0.0);
  alphatot = 1.0/alphatot;
  double alpha[4];
  for (int i=0; i<4; i++){
    alpha[i] = max(0.0, dphii[i]/dphi)*alphatot;
    //fprintf(stdout, "alpha[ilocal=%d] = %e\n", i, alpha[i]);
  }
  if(debug){
    fprintf(stdout, "dphi = %e - alpha = %e %e %e %e\n", dphi, alpha[0],alpha[1],alpha[2],alpha[3]);
  }


  //update PsiRes and w and beta
  for (int i=0; i<4; i++){
    locdphi[i] += alpha[i]*(dphi - Fsrc*vol);
    locw[i] += alpha[i]*vol;
    if (alpha[i]>0.0)
      locbeta[i] += Kp[i]/alphatot;
  }

}
//------------------------------------------------------------------------------
double Tet::findRootPolynomialNewtonRaphson(double f1, double f2, double fp1, double fp2)
{
// finds one root between 0 and 1 of the Hermite interpolation polynomial
// that verifies P(0)  = f1    P(1)  = f2
//               P'(0) = fp1   P'(1) = fp2

  double coeff[4] = { 2.0*(f1-f2)+fp1+fp2, -3.0*(f1-f2)-2.0*fp1-fp2, fp1, f1};
  double coeffp[3] = {3.0*coeff[0],2.0*coeff[1],coeff[2]};
  
  double eps = 1.e-6;                  //precision
  double xn = 0.5;                     //initial guess
  bool notConverged = true;
  int maxIts = 100;
  int it = 0;
  int ierr = 0;
  double xnp1, f, fp, xn2,xn3;
  while(notConverged){
    xn2 = xn*xn;
    xn3 = xn*xn2;
    f  = coeff[0]*xn3 + coeff[1]*xn2 + coeff[2]*xn + coeff[3];
    fp = coeffp[0]*xn2 + coeffp[1]*xn + coeffp[2];
    assert(fp!= 0.0);
    xnp1 = xn - f/fp;
    if( fabs((xnp1-xn)/(xnp1+xn)) < eps) notConverged = false;
    xn = xnp1;
    it++;
    if(it>maxIts){
      ierr++;
      fprintf(stdout, "*** Error: max iteration reached in Newton-Raphson solver for Hermite\n");
      notConverged = false;
    }
  }

  // check value of xn
  if(xn<0.0 || xn>1.0) {
    ierr++;
    fprintf(stdout, "*** Error: solution(%e) is out of bound in Hermite polynomial root finder\n",xn);
  }

  return xn;
}
//------------------------------------------------------------------------------
extern int zroots(bcomp *a, int degree, bcomp *roots, const bool &polish);
int Tet::findRootPolynomialLaguerre(double f1, double f2, double fp1, double fp2,
                                    double &root)
{
/* The Laguerre method (cf Numerical Recipes in C++) is used to find
** the roots of the polynomial we are considering.
** However we need only one solution that lies between 0 and 1.
** We choose the one that fits, and if there are several of them, 
** we revert to a linear interpolation to find the interface location
** instead of Hermite interpolation polynomial.
** Function returns the number of real roots in [0,1]. If there are
** several, root contains the 'last' corresponding one.
*/

  int degree = 3; //degree of the polynomial
  bcomp coeff[degree+1]; //coeff of the polynomial
  coeff[3] = bcomp(2.0*(f1-f2)+fp1+fp2,0.0);
  coeff[2] = bcomp(-3.0*(f1-f2)-2.0*fp1-fp2,0.0);
  coeff[1] = bcomp(fp1,0.0);
  coeff[0] = bcomp(f1,0.0);


  bcomp roots[degree]; //roots of the polynomial
  int err = zroots(coeff, degree, roots, true);
  if(err>0) return 1000;

  //check which real one is in the bounds [0,1]
  int counter = 0;
  int index = -1;
  double eps = 1.0e-7;
  for(int i=0; i<degree; i++)
    if(fabs(imag(roots[i]))<eps*fabs(real(roots[i])) &&  // check that the root is real
       real(roots[i]) <= 1.0                         &&  // check bounds of that root
       real(roots[i]) >= 0.0                          ){
      counter++;
      index = i;
    }

  root = real(roots[index]); //we assume f1*f2<0 and 
                             //thus there must be such a solution 
                             //and index should be well defined!
  return counter;
}

//------------------------------------------------------------------------------
bool Tet::computeDistancePlusPhiToOppFace(double phi[3], Vec3D Y0,
                                          Vec3D Y1, Vec3D Y2, double &mini, bool show)
{
  bool found  = false;
  double eps = 1.0e-14;
  double one  = 1.0;
  double zero = 0.0;
	double tol1 = eps*Y1.norm();
	double tol2 = eps*Y2.norm();
  // mini is assumed to already have some distance, usually set by 
  // distance to vertices
	
  // strategy is as follows:
  // find point on opposite face where minimum is attained
  // the point to find is defined by a vector originating on the node i
  // this vector Z is equal to Z_ortho+Z2*y2+Z3*y3
  // where Z_ortho is the orthogonal component wrt to the opp face
  //       y2 is the vector Y1 (edge between nodes 0 and 1 of opp face)
  //       y3 is the vector Y2 (edge between nodes 0 and 2 of opp face)
  // This vector can be expressed as Y0 + alpha*Y1 + beta*Y2
  // We solve a 2-by-2 system to get alpha and beta.
  // Finally, we compute values at that point and thus determine the minimum


  Vec3D normalY = Y1 ^ Y2;
  normalY /= normalY.norm();
  double orthogonal = Y0 * normalY;

  // constants
  double y2sq = Y1 * Y1;
  double y3sq = Y2 * Y2;
  double K  = Y1 * Y2;
  // to solve 2-by-2 system
  double det = y2sq*y3sq - K*K;
  det = 1.0/det;
  double alpha, beta;

  if(fabs(phi[1])<tol1 && fabs(phi[2])<tol2){
    //solution is the projection of X[nodeNum[i]] on that face
    double rhs[2] = {Y0 * Y1, Y0 * Y2};
    alpha = det *(y3sq*rhs[0]-   K*rhs[1]);
    beta  = det *(  -K*rhs[0]+y2sq*rhs[1]);
    if(alpha>= zero && alpha <= one &&
       beta >= zero && beta  <= one &&
       alpha+beta<= one){
			if(show) fprintf(stdout, "face1 - %e\n", phi[0]+alpha*phi[1]+beta*phi[2] + (Y0 - alpha*Y1 - beta*Y2).norm()); 
      mini = min(mini,phi[0]+alpha*phi[1]+beta*phi[2] + (Y0 - alpha*Y1 - beta*Y2).norm());
      return true;
    }
    return false;
  }else if(fabs(phi[1])<tol1){
    //
    double temp = y3sq - K*K/y2sq;
    //if(fabs(temp-phi[2]*phi[2])<tol2*tol2) return false;
		assert(temp-phi[2]*phi[2]!=0.0);
    double Z3sq = phi[2]*phi[2]*orthogonal*orthogonal/(temp*(temp-phi[2]*phi[2]));
    if(Z3sq<0.0) return false;

    //first solution
    double Z3 = sqrt(Z3sq);
    double Z2 = -K*Z3/y2sq;
    double rhs[2] = {Y0*Y1 - (Z2*y2sq+Z3*K), Y0*Y2 - (Z3*y3sq+Z2*K)};
    alpha = det *(y3sq*rhs[0]-   K*rhs[1]);
    beta  = det *(  -K*rhs[0]+y2sq*rhs[1]);
    if(show) fprintf(stdout, "case1       alpha = %e and beta = %e\n", alpha, beta);
    if(alpha>= zero && alpha <= one &&
       beta >= zero && beta  <= one &&
       alpha+beta<= one){
			if(show) fprintf(stdout, "face2 - %e\n", phi[0]+alpha*phi[1]+beta*phi[2] + (Y0 - alpha*Y1 - beta*Y2).norm()); 
      mini = min(mini,phi[0]+alpha*phi[1]+beta*phi[2] + (Y0 - alpha*Y1 - beta*Y2).norm());
      found = true;
    }
    //second solution
    Z3 = -Z3;
    Z2 = -K*Z3/y2sq;
    rhs[0] = Y0*Y1 - (Z2*y2sq+Z3*K);
    rhs[1] = Y0*Y2 - (Z3*y3sq+Z2*K);
    alpha = det *(y3sq*rhs[0]-   K*rhs[1]);
    beta  = det *(  -K*rhs[0]+y2sq*rhs[1]);
    if(show) fprintf(stdout, "case1       alpha = %e and beta = %e\n", alpha, beta);
    if(alpha>= zero && alpha <= one &&
       beta >= zero && beta  <= one &&
       alpha+beta<= one){
      mini = min(mini,phi[0]+alpha*phi[1]+beta*phi[2] + (Y0 - alpha*Y1 - beta*Y2).norm());
			if(show) fprintf(stdout, "face3 - %e\n", phi[0]+alpha*phi[1]+beta*phi[2] + (Y0 - alpha*Y1 - beta*Y2).norm()); 
      found = true;
    }
    return found;

  }else if(fabs(phi[2])<tol2){
    //
    double temp = y2sq - K*K/y3sq;
    //if(fabs(temp-phi[1]*phi[1])<tol1*tol1) return false;
		assert(temp-phi[1]*phi[1]!=0.0);
    double Z2sq = phi[1]*phi[1]*orthogonal*orthogonal/(temp*(temp-phi[1]*phi[1]));
    if(Z2sq<0.0) return false;

    //first solution
    double Z2 = sqrt(Z2sq);
    double Z3 = -K*Z2/y3sq;
    double rhs[2] = {Y0*Y1 - (Z2*y2sq+Z3*K), Y0*Y2 - (Z3*y3sq+Z2*K)};
    alpha = det *(y3sq*rhs[0]-   K*rhs[1]);
    beta  = det *(  -K*rhs[0]+y2sq*rhs[1]);
    if(show) fprintf(stdout, "case2       alpha = %e and beta = %e\n", alpha, beta);
    if(alpha>= zero && alpha <= one &&
       beta >= zero && beta  <= one &&
       alpha+beta<= one){
      mini = min(mini,phi[0]+alpha*phi[1]+beta*phi[2] + (Y0 - alpha*Y1 - beta*Y2).norm());
			if(show) fprintf(stdout, "face4 - %e\n", phi[0]+alpha*phi[1]+beta*phi[2] + (Y0 - alpha*Y1 - beta*Y2).norm()); 
      found = true;
    }
    //second solution
    Z2 = -Z2;
    Z3 = -K*Z2/y3sq;
    rhs[0] = Y0*Y1 - (Z2*y2sq+Z3*K);
    rhs[1] = Y0*Y2 - (Z3*y3sq+Z2*K);
    alpha = det *(y3sq*rhs[0]-   K*rhs[1]);
    beta  = det *(  -K*rhs[0]+y2sq*rhs[1]);
    if(show) fprintf(stdout, "case2       alpha = %e and beta = %e\n", alpha, beta);
    if(alpha>= zero && alpha <= one &&
       beta >= zero && beta  <= one &&
       alpha+beta<= one){
			if(show) fprintf(stdout, "face5 - %e\n", phi[0]+alpha*phi[1]+beta*phi[2] + (Y0 - alpha*Y1 - beta*Y2).norm()); 
      mini = min(mini,phi[0]+alpha*phi[1]+beta*phi[2] + (Y0 - alpha*Y1 - beta*Y2).norm());
      found = true;
    }
    return found;

  }else{
    //general case
    double K2 = phi[1]*K - phi[2]*y2sq;
    double K3 = phi[2]*K - phi[1]*y3sq;
    double K2sq = K2*K2;
    double K3sq = K3*K3;
    double denom = -phi[1]*phi[2]*(K3sq*y2sq+K2sq*y3sq+2.0*K*K2*K3);
    denom += (K*K3+y3sq*K2)*(K3*y2sq+K*K2);
		if(show) fprintf(stdout, "tol1 = %e -- tol2 = %e -- K2 = %e -- K3 = %e -- denom = %e\n", tol1,tol2,K2, K3, denom);
		//if(fabs(denom)<tol1*tol2) return false;
		if(denom==0.0) return false;

    double Z2,Z3;
    //first solution
    if(K3!=0.0){
      double Z2sq = (K3sq*phi[1]*phi[2]*orthogonal*orthogonal)/denom;
      if(Z2sq<0.0) return false;
      Z2 = sqrt(Z2sq);
      Z3 = K2*Z2/K3;
    }else if(K2!=0.0){
      double Z3sq = (K2sq*phi[1]*phi[2]*orthogonal*orthogonal)/denom;
      if(Z3sq<0.0) return false;
      Z3 = sqrt(Z3sq);
      Z2 = K3*Z3/K2;
    }else{
      fprintf(stdout, "***Error: K3 = K2 = 0 in reinitialization\n");
      exit(1);
    }
    double rhs[2] = {Y0*Y1 - (Z2*y2sq+Z3*K), Y0*Y2 - (Z3*y3sq+Z2*K)};
    alpha = det *(y3sq*rhs[0]-   K*rhs[1]);
    beta  = det *(  -K*rhs[0]+y2sq*rhs[1]);
    if(show) fprintf(stdout, "case3       alpha = %e and beta = %e\n", alpha, beta);
    if(alpha>= zero && alpha <= one &&
       beta >= zero && beta  <= one &&
       alpha+beta<= one){
			if(show) fprintf(stdout, "face6 - %e %e %e\n", alpha, beta, phi[0]+alpha*phi[1]+beta*phi[2] + (Y0 - alpha*Y1 - beta*Y2).norm()); 
      mini = min(mini,phi[0]+alpha*phi[1]+beta*phi[2] + (Y0 - alpha*Y1 - beta*Y2).norm());
      found = true;
    }

    //second solution
    Z2 = -Z2;
    Z3 = -Z3;
    rhs[0] = Y0*Y1 - (Z2*y2sq+Z3*K);
    rhs[1] = Y0*Y2 - (Z3*y3sq+Z2*K);
    alpha = det *(y3sq*rhs[0]-   K*rhs[1]);
    beta  = det *(  -K*rhs[0]+y2sq*rhs[1]);
    if(show) fprintf(stdout, "case3       alpha = %e and beta = %e\n", alpha, beta);
    if(alpha>= zero && alpha <= one &&
       beta >= zero && beta  <= one &&
       alpha+beta<= one){
			if(show) fprintf(stdout, "face7 - %e\n", phi[0]+alpha*phi[1]+beta*phi[2] + (Y0 - alpha*Y1 - beta*Y2).norm()); 
      mini = min(mini,phi[0]+alpha*phi[1]+beta*phi[2] + (Y0 - alpha*Y1 - beta*Y2).norm());
      found = true;
    }
    return found;
  }

  return false;

}
//------------------------------------------------------------------------------
bool Tet::computeDistancePlusPhiToEdges(double phi[3], Vec3D Y0,
                                        Vec3D Y1, Vec3D Y2, double &mini, bool show)
{
  bool found = false;
  //1st edge
  found = computeDistancePlusPhiToEdge(phi[0],phi[1],Y0,Y1,mini, show);

  //2nd edge
  found = computeDistancePlusPhiToEdge(phi[0],phi[2],Y0,Y2,mini, show);

  //3rd edge
  found = computeDistancePlusPhiToEdge(phi[1]+phi[0],phi[2]-phi[1],Y0-Y1,Y2-Y1,mini, show);

  return found;
}
//------------------------------------------------------------------------------
bool Tet::computeDistancePlusPhiToEdge(double phi0, double phi1,
                                       Vec3D Y0, Vec3D Y1, double &mini, bool show)
{
  bool found = false;
	double eps = 1.0e-14;
  // same approach is used as in computeDistancePlusPhiToOppFace
  // except that it is one dimensional (along an edge) and thus we 
  // have only one parameter alpha (instead of two).

  // normal to Y1 in (Y0,Y1)-plane can be written n = n0*Y0+n1*Y1 
  // note that Y0 and Y1 are not necessarily orthogonal.


  double K = Y0 * Y1;
  double y1sq = Y1 * Y1;
  Vec3D n = Y0 -K/y1sq * Y1; 
	double normn = n.norm();
	if(normn==0.0) return false;
  n /= normn;
  double orthogonal = Y0 * n;

  if(y1sq==phi1*phi1) return false;
	//if(fabs(y1sq-phi1*phi1)<eps*Y1.norm()) return false;
  double Z1sq = phi1*phi1*orthogonal*orthogonal/(y1sq*(y1sq-phi1*phi1));
  if(Z1sq<0.0) return false;

  double Z1 = sqrt(Z1sq);
  double alpha = (K - Z1*y1sq)/y1sq;
  if(alpha>=0.0 && alpha<=1.0){
	  if(show) fprintf(stdout, "edge1 -%e %e\n", alpha, phi0+alpha*phi1+ (Y0 - alpha*Y1).norm()); 
    found = true;
    mini = min(mini,phi0+alpha*phi1 + (Y0 - alpha*Y1).norm());
  }
  Z1 = -Z1;
  alpha = (K - Z1*y1sq)/y1sq;
  if(alpha>=0.0 && alpha<=1.0){
    found = true;
	  if(show) fprintf(stdout, "edge2 -%e %e\n", alpha, phi0+alpha*phi1+ (Y0 - alpha*Y1).norm()); 
    mini = min(mini,phi0+alpha*phi1 + (Y0 - alpha*Y1).norm());
  }

  return found;
}
//------------------------------------------------------------------------------
bool Tet::computeDistancePlusPhiToVertices(double phi[3], Vec3D Y0,
                                           Vec3D Y1, Vec3D Y2, double &mini, bool show)
{
  // only possibilities left are the vertices.
  mini = phi[0]+Y0.norm();
  mini = min(mini, phi[1]+phi[0]+(Y0-Y1).norm());
  mini = min(mini, phi[2]+phi[0]+(Y0-Y2).norm());
	if(show) fprintf(stdout, "vertices - %e %e %e\n", phi[0]+Y0.norm(),phi[1]+phi[0]+(Y0-Y1).norm(),phi[2]+phi[0]+(Y0-Y2).norm());
	//if(mini<1.0e-14*Y0.norm() || mini<1.0e-14*Y1.norm() || mini<1.0e-14*Y2.norm()) return true;
	return false;

}
//------------------------------------------------------------------------------
int Tet::computeDistanceToAll(double phi[3], Vec3D Y0, Vec3D Y1, Vec3D Y2, double &psi)
{
  double eps = 1.0e-14;

  // psi is overwritten here by this function!
  computeDistancePlusPhiToVertices(phi,Y0,Y1,Y2,psi);


  bool found = computeDistancePlusPhiToOppFace(phi,Y0,Y1,Y2,psi);
	if(found) return 1;

  computeDistancePlusPhiToEdges(phi,Y0,Y1,Y2,psi);
	//if(psi<eps*Y0.norm()||psi<eps*Y1.norm()||psi<eps*Y2.norm()) return 1;
  return -1;

}


