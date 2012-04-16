#ifndef _POLYGONRECONSTRUCTIONDATA_H_
#define _POLYGONRECONSTRUCTIONDATA_H_

class Elem;
class LevelSetStructure;
template<class Scalar, int dim> class SVec;

struct Vec3D;

//------------------------------------------------------------------------------

struct PolygonReconstructionData { //for force computation under the embedded framework
    PolygonReconstructionData() : numberOfEdges(0) {}
    int numberOfEdges;
    int edgeWithVertex[4][2];
    int edge[4];

    void AssignSingleEdge(int n1, int n2,int l1){ //for PhysBAM only
        numberOfEdges=1;
        edge[0]=l1;
        edgeWithVertex[0][0]=n1; edgeWithVertex[0][1]=n2;
    }

    void AssignTwoEdges(int n1, int n2, int n3, int l1, int l2){ //for PhysBAM only
        numberOfEdges=2;
        edge[0]=l1; edge[1]=l2;
        edgeWithVertex[0][0]=n1; edgeWithVertex[0][1]=n2;
        edgeWithVertex[1][0]=n1; edgeWithVertex[1][1]=n3;
    }

    void AssignTriangle(int n1, int n2, int n3, int n4, int l1, int l2, int l3, bool owned_by_single_vertex=true){
        numberOfEdges=3;
        edge[0]=l1; edge[1]=l2; edge[2]=l3;
        if(owned_by_single_vertex){
            edgeWithVertex[0][0]=n1; edgeWithVertex[0][1]=n2;
            edgeWithVertex[1][0]=n1; edgeWithVertex[1][1]=n3;
            edgeWithVertex[2][0]=n1; edgeWithVertex[2][1]=n4;}
        else{
            edgeWithVertex[0][0]=n2; edgeWithVertex[0][1]=n1;
            edgeWithVertex[1][0]=n3; edgeWithVertex[1][1]=n1;
            edgeWithVertex[2][0]=n4; edgeWithVertex[2][1]=n1;}
    }

    void AssignQuadTriangle(int n1, int n2, int n3, int n4, int l1, int l2, int l3){ //for PhysBAM only
        numberOfEdges=3;
        edge[0]=l1; edge[1]=l2; edge[2]=l3;
        edgeWithVertex[0][0]=n1; edgeWithVertex[0][1]=n3;
        edgeWithVertex[1][0]=n1; edgeWithVertex[1][1]=n4;
        edgeWithVertex[2][0]=n2; edgeWithVertex[2][1]=n3;
    }

    void AssignQuadrilateral(int n1, int n2, int n3, int n4, int l1, int l2, int l3, int l4){
        numberOfEdges=4;
        edge[0]=l1; edge[1]=l2; edge[2]=l3; edge[3]=l4;
        edgeWithVertex[0][0]=n1; edgeWithVertex[0][1]=n3;
        edgeWithVertex[1][0]=n1; edgeWithVertex[1][1]=n4;
        edgeWithVertex[2][0]=n2; edgeWithVertex[2][1]=n4;
        edgeWithVertex[3][0]=n2; edgeWithVertex[3][1]=n3;
    }
};

//------------------------------------------------------------------------------

int getPolygons(Elem &elem, LevelSetStructure &LSS, PolygonReconstructionData* polygons);
void getPolygonNormal(SVec<double,3>& X, Vec3D &normal, LevelSetStructure &LSS, PolygonReconstructionData &polygon);

#endif
