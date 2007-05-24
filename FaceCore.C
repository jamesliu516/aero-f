#include <Face.h>

#include <RefVal.h>
#include <BcDef.h>
#include <Edge.h>
#include <Vector3D.h>
#include <Vector.h>
#include <BinFileHandler.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef OLD_STL
#include <defalloc.h>
#include <algo.h>
#else
#include <algorithm>
using std::stable_sort;
using std::min;
using std::max;
using std::swap;
#endif

const double Face::third = 1.0/3.0;
const int Face::edgeEnd[3][2] = { {0,1}, {1,2}, {2,0} };

//------------------------------------------------------------------------------

FaceSet::FaceSet(int value)  
{

  numFaces = value;
  faces	= new Face[value];

}

//------------------------------------------------------------------------------

void Face::setup(int fc, int *nn, int sid)
{

  if (fc == BC_OUTLET_MOVING || fc == BC_OUTLET_FIXED || 
      fc == BC_INLET_MOVING || fc == BC_INLET_FIXED ||
      fc == BC_ADIABATIC_WALL_MOVING || fc == BC_ADIABATIC_WALL_FIXED ||
      fc == BC_SLIP_WALL_MOVING || fc == BC_SLIP_WALL_FIXED ||
      fc == BC_ISOTHERMAL_WALL_MOVING || fc == BC_ISOTHERMAL_WALL_FIXED ||
      fc == BC_SYMMETRY)
    code = fc;
  else {
    fprintf(stderr, "*** Error: incorrect boundary face code (%d)\n", fc);
    exit(1);
  }

  nodeNum[0] = nn[0];
  nodeNum[1] = nn[1];
  nodeNum[2] = nn[2];

  surface_id = sid;
}

//------------------------------------------------------------------------------

void Face::setType(int *facemap)
{

  code = facemap[code];

}

//------------------------------------------------------------------------------

void Face::setNodeType(int* priority, int* nodeType)
{

  for (int j=0; j<3; ++j)
    if (priority[code] > priority[ nodeType[ nodeNum[j] ] ])
      nodeType[ nodeNum[j] ] = code;

}

//------------------------------------------------------------------------------

void Face::setNodeFaceType(int* nodeFaceType)
{

  int type;

  if(code==BC_INLET_FIXED  || code==BC_INLET_MOVING ||
     code==BC_OUTLET_FIXED || code==BC_OUTLET_MOVING ){
    for (int j=0; j<3; j++){
      type = nodeFaceType[nodeNum[j]];
      if(type == 0)
        nodeFaceType[nodeNum[j]] = 1;
      if(type == -1)
        nodeFaceType[nodeNum[j]] = 2;
    }
  }else{ //only possibilities left should be wall and symmetry
    for (int j=0; j<3; j++){
      type = nodeFaceType[nodeNum[j]];
      if(type == 0)
        nodeFaceType[nodeNum[j]] = -1;
      if(type == 1)
        nodeFaceType[nodeNum[j]] = 2;
    }
  }

}

//------------------------------------------------------------------------------

void Face::setElementNumber(int num, int* nn)
{

  elemNum = num;
  if (nn) {
    nodeNum[0] = nn[0];
    nodeNum[1] = nn[1];
    nodeNum[2] = nn[2];
  }

}

//------------------------------------------------------------------------------

void Face::tagNodesOnBoundaries(Vec<bool> &tagNodes)
{

  tagNodes[ nodeNum[0] ] = true;
  tagNodes[ nodeNum[1] ] = true;
  tagNodes[ nodeNum[2] ] = true;

}

//------------------------------------------------------------------------------

void Face::tagEdgesOnBoundaries(Vec<bool> &tagEdges)
{

  tagEdges[ edgeNum[0] ] = true;
  tagEdges[ edgeNum[1] ] = true;
  tagEdges[ edgeNum[2] ] = true;

}

//------------------------------------------------------------------------------

void Face::reorder()
{ 

#ifdef OLD_STL
  sort(&nodeNum[0], &nodeNum[3]); 
#else
  stable_sort(&nodeNum[0], &nodeNum[3]); 
#endif

}

//------------------------------------------------------------------------------

void Face::numberEdges(EdgeSet &edges)
{

  int numEdges = edges.size();

  edgeNum[0] = edges.find(nodeNum[ edgeEnd[0][0] ], nodeNum[ edgeEnd[0][1] ]);
  edgeNum[1] = edges.find(nodeNum[ edgeEnd[1][0] ], nodeNum[ edgeEnd[1][1] ]);
  edgeNum[2] = edges.find(nodeNum[ edgeEnd[2][0] ], nodeNum[ edgeEnd[2][1] ]);

  if (edges.size() != numEdges) {
    fprintf(stderr, "*** Error: could not find face egdes\n");
    exit(1);
  }

}

//------------------------------------------------------------------------------

double Face::computeVolume(Vec3D &xa_n, Vec3D &xb_n, Vec3D &xc_n, 
			   Vec3D &xa_np1, Vec3D &xb_np1, Vec3D &xc_np1)
{

  const static double sixth = 1.0/6.0;

  Vec3D dxa = xa_np1 - xa_n;
  Vec3D dxb = xb_np1 - xb_n;
  Vec3D dxc = xc_np1 - xc_n;
  Vec3D xac_n = xa_n - xc_n;
  Vec3D xbc_n = xb_n - xc_n;
  Vec3D xac_np1 = xa_np1 - xc_np1;
  Vec3D xbc_np1 = xb_np1 - xc_np1;
  
  Vec3D eta = sixth * ((xac_n ^ xbc_n) + (xac_np1 ^ xbc_np1) + 
		       0.5 * ((xac_n ^ xbc_np1) + (xac_np1 ^ xbc_n)));
  
  return third * (dxa + dxb + dxc) * eta;

}

//------------------------------------------------------------------------------
// computation of the OUTWARD face normal

void Face::computeNormal(SVec<double,3> &X, Vec3D &faceNorm)
{

  Vec3D x[3] = {X[ nodeNum[0] ], X[ nodeNum[1] ], X[ nodeNum[2] ]};

  faceNorm = 0.5 * ((x[2] - x[0]) ^ (x[1] - x[0]));

}

//------------------------------------------------------------------------------

void Face::computeEdgeNormals(SVec<double,3>& X, int* l2gl, SVec<double,6>& normals)
{

  Vec3D faceNorm;
  computeNormal(X, faceNorm);
  faceNorm *= -1.0 / sqrt(faceNorm*faceNorm);

  for (int l=0; l<3; ++l) {
    bool ori;
    if (l2gl[ nodeNum[ edgeEnd[l][0] ] ] < l2gl[ nodeNum[ edgeEnd[l][1] ] ])
      ori = true;
    else
      ori = false;

    if (ori) {
      normals[ edgeNum[l] ][0] = faceNorm[0];
      normals[ edgeNum[l] ][1] = faceNorm[1];
      normals[ edgeNum[l] ][2] = faceNorm[2];
    }
    else {
      normals[ edgeNum[l] ][3] = faceNorm[0];
      normals[ edgeNum[l] ][4] = faceNorm[1];
      normals[ edgeNum[l] ][5] = faceNorm[2];
    }
  }

}

//------------------------------------------------------------------------------
// computation of the OUTWARD face normal

void Face::computeNormalGCL1(SVec<double,3> &Xn, SVec<double,3> &Xnp1, 
			     SVec<double,3> &Xdot, Vec3D &faceNorm, 
			     double &faceNormVel)
{

  static double twelfth = 1.0/12.0;

  Vec3D x_n[3] = {Xn[ nodeNum[0] ], Xn[ nodeNum[1] ], Xn[ nodeNum[2] ]};
  Vec3D x_np1[3] = {Xnp1[ nodeNum[0] ], Xnp1[ nodeNum[1] ], Xnp1[ nodeNum[2] ]};
  Vec3D xdot[3] = {Xdot[ nodeNum[0] ], Xdot[ nodeNum[1] ], Xdot[ nodeNum[2] ]};

  Vec3D x01_n = x_n[1] - x_n[0];
  Vec3D x02_n = x_n[2] - x_n[0];

  Vec3D x01_np1 = x_np1[1] - x_np1[0];
  Vec3D x02_np1 = x_np1[2] - x_np1[0];

  faceNorm = twelfth * (((2.0*x02_np1 + x02_n) ^ x01_np1) + 
			((2.0*x02_n + x02_np1) ^ x01_n));

  faceNormVel = third * (xdot[0] + xdot[1] + xdot[2]) * faceNorm;

}

//------------------------------------------------------------------------------
// computation of the OUTWARD face normal

void Face::computeNormalEZGCL1(double oodt, SVec<double,3> &Xn, SVec<double,3> &Xnp1, 
			       Vec3D &faceNorm, double &faceNormVel)
{

  Vec3D x_n[3] = {Xn[ nodeNum[0] ], Xn[ nodeNum[1] ], Xn[ nodeNum[2] ]};
  Vec3D x_np1[3] = {Xnp1[ nodeNum[0] ], Xnp1[ nodeNum[1] ], Xnp1[ nodeNum[2] ]};

  faceNorm = 0.5 * ((x_np1[2] - x_np1[0]) ^ (x_np1[1] - x_np1[0]));
  //EZ1 faceNorm = 0.5 * ((x_n[2] - x_n[0]) ^ (x_n[1] - x_n[0]));

  double vol = computeVolume(x_n[0], x_n[2], x_n[1], x_np1[0], x_np1[2], x_np1[1]);

  faceNormVel = oodt * vol;

}

//------------------------------------------------------------------------------

int FaceSet::read(BinFileHandler &file, int numRanges, int (*ranges)[2], int *map)
{
 
  // read in number of faces in cluster (not used)
  int numClusFaces;
  file.read(&numClusFaces, 1);

  // read in the offset for the first face
  BinFileHandler::OffType start, tocStart;
  file.read(&start, 1);

  // set the location of first face in the TOC
  tocStart = file.tell();

  int count = 0;

  // read in ranges
  int nSymm = 0;
  for (int iRange = 0; iRange < numRanges; ++iRange) {

    // compute number of faces in range
    int nFaces = ranges[iRange][1] - ranges[iRange][0] + 1;

    // seek to correct position in file thanks to the table of contents
    file.seek(tocStart + sizeof(int) * ranges[iRange][0]);
    int toc;
    file.read(&toc, 1);
    file.seek(start + sizeof(int) * toc);
    for (int i = 0; i < nFaces; i++)  {
      // read in the face type (not used)
      int type;
      file.read(&type, 1);
      // read in global face number
      file.read(map + count, 1);
      // read in face nodes
      int nodeNum[3];
      file.read(nodeNum, 3);
      // read in the face code
      int code;
      file.read(&code, 1);
      if (code == 6)  nSymm++;
      // PJSA read in the surface id
      int surface_id;
      file.read(&surface_id, 1);
      faces[count].setup(code, nodeNum, surface_id);
      count++;
    }
  }

  if (count != numFaces) {
    fprintf(stderr, "*** Error: wrong number of faces read (%d instead of %d)\n",
	    count, numFaces);
    exit(1);
  }

  return numClusFaces;

}

