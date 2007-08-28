#include <Elem.h>

#include <Edge.h>
#include <Face.h>
#include <Vector3D.h>
#include <Vector.h>
#include <BinFileHandler.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//------------------------------------------------------------------------------

double Elem::computeLongestEdge(SVec<double,3> &X)
{
  double result = 0.0;
  double temp;

  for (int i=0; i<numEdges(); i++) {
    Vec3D x0 = X[nodeNum(edgeEnd(i,0))];
    Vec3D x1 = X[nodeNum(edgeEnd(i,1))];
    Vec3D dx = x0-x1;

    double temp = dx.norm();

    if (temp > result) result = temp;      
  }

  return result;
}

//------------------------------------------------------------------------------

void Elem::numberEdges(EdgeSet &edges)
{
  for (int i=0; i<numEdges(); i++)
    edgeNum(i) = edges.find(nodeNum( edgeEnd(i,0) ), nodeNum( edgeEnd(i,1) ));
}

//------------------------------------------------------------------------------

void Elem::renumberEdges(Vec<int> &newNum)
{
  for (int l=0; l<numEdges(); ++l) 
    edgeNum(l) = newNum[ edgeNum(l) ];
}

//------------------------------------------------------------------------------

int Elem::countNodesOnBoundaries(Vec<bool> &tagNodes)
{
  int num = 0;

  for (int i=0; i<numNodes(); i++)
    if (tagNodes[ nodeNum(i) ]) num++;

  return num;
}

//------------------------------------------------------------------------------

int Elem::setFaceToElementConnectivity(int i, Vec<bool> &tagNodes, 
					  MapFaces &mf, int(*fm)[2])
{
  int nswap = 0;
  
  for (int k=0; k<numFaces(); k++) {
    int nn[Face::MaxNumNd];
    
    for (int l=0; l<faceNnd(k); l++)
      nn[l] = nodeNum( faceDef(k,l) );
    
    // This is not the complete test, but it's no big deal if we check
    // for more faces than necessary.
    if (tagNodes[ nn[0] ] && tagNodes[ nn[1] ] && tagNodes[ nn[2] ]) {
      MaxFace f(faceNnd(k), nn);
      f.reorder();
      
      MapFaces::iterator it = mf.find(f);
      if (it != mf.end()) {
	int fn = (*it).second;
	fm[fn][0] = i;
	fm[fn][1] = f.rotDir * it->first.rotDir;
	if(fm[fn][1] < 0)
	  nswap++;
      }
    }
  }
  
  return nswap;
}

//------------------------------------------------------------------------------

ElemSet::ElemSet(int value)  
{

  numElems = value;
  elems = new Elem*[value];

}

//------------------------------------------------------------------------------

ElemSet::~ElemSet()  
{
  
  if (elems) delete[] elems;

}

//------------------------------------------------------------------------------

int ElemSet::read(BinFileHandler &file, int numRanges, int (*ranges)[2], int *map)  
{

  // read in number of elems in cluster (not used)
  int numClusElems;
  file.read(&numClusElems, 1);

  // read in the offset for the first elem
  BinFileHandler::OffType start, tocStart;
  file.read(&start, 1);

  // set the location of first elem in the TOC
  tocStart = file.tell();

  int count = 0;

  // read in ranges
  for (int iRange = 0; iRange < numRanges; ++iRange) {

    // compute number of elems in range
    int nElems = ranges[iRange][1] - ranges[iRange][0] + 1;

    // seek to correct position in file thanks to the table of contents
    file.seek(tocStart + sizeof(int) * ranges[iRange][0]);
    int toc;
    file.read(&toc, 1);
    file.seek(start + sizeof(int) * toc);

    for (int i = 0; i < nElems; i++)  {
      int type, volume_id;

      // read in the element type
      file.read(&type, 1);

      // create elem in elemset depending on type
      switch (type) {
      case Elem::TET:
	addElem(count, new(memElems) ElemTet);
	break;

      default:
	fprintf(stderr, "Error: Could not find elem type %d\n", type);
	exit(1);
      }

      // read in global elem number
      file.read(&map[count], 1);

      // read in volume id (for porous)
      file.read(&volume_id, 1);
      elems[count]->setVolumeID(volume_id);

      // read in elem nodes
      file.read( elems[count]->nodeNum(), elems[count]->numNodes());

      // count number of elems
      count++;
    }
  }

  if (count != numElems) {
    fprintf(stderr, "*** Error: wrong number of elems read (%d instead of %d)\n",
	    count, numElems);
    exit(1);
  }

  return numClusElems;

}

//------------------------------------------------------------------------------

void ElemSet::computeDynamicLESTerm(DynamicLESTerm *dles, SVec<double,2> &CsDeltaSq, SVec<double,3> &X,
				    Vec<double> &Cs, Vec<double> &VolSum)
{

 for (int i=0; i<numElems; ++i)
    elems[i]->computeDynamicLESTerm(dles, CsDeltaSq, X, Cs, VolSum);

}

//------------------------------------------------------------------------------