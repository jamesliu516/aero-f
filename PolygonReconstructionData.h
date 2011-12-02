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

    void AssignSingleEdge(int n1, int n2){ //for PhysBAM only
        numberOfEdges=1;
        edgeWithVertex[0][0]=n1; edgeWithVertex[0][1]=n2;
    }

    void AssignTwoEdges(int n1, int n2, int n3){ //for PhysBAM only
        numberOfEdges=2;
        edgeWithVertex[0][0]=n1; edgeWithVertex[0][1]=n2;
        edgeWithVertex[1][0]=n1; edgeWithVertex[1][1]=n3;
    }

    void AssignTriangle(int n1, int n2, int n3, int n4,bool owned_by_single_vertex=true){
        numberOfEdges=3;
        if(owned_by_single_vertex){
            edgeWithVertex[0][0]=n1; edgeWithVertex[0][1]=n2;
            edgeWithVertex[1][0]=n1; edgeWithVertex[1][1]=n3;
            edgeWithVertex[2][0]=n1; edgeWithVertex[2][1]=n4;}
        else{
            edgeWithVertex[0][0]=n2; edgeWithVertex[0][1]=n1;
            edgeWithVertex[1][0]=n3; edgeWithVertex[1][1]=n1;
            edgeWithVertex[2][0]=n4; edgeWithVertex[2][1]=n1;}
    }

    void AssignQuadTriangle(int n1, int n2, int n3, int n4){ //for PhysBAM only
        numberOfEdges=3;
        edgeWithVertex[0][0]=n1; edgeWithVertex[0][1]=n3;
        edgeWithVertex[1][0]=n1; edgeWithVertex[1][1]=n4;
        edgeWithVertex[2][0]=n2; edgeWithVertex[2][1]=n3;
    }

    void AssignQuadrilateral(int n1, int n2, int n3, int n4){
        numberOfEdges=4;
        edgeWithVertex[0][0]=n1; edgeWithVertex[0][1]=n3;
        edgeWithVertex[1][0]=n1; edgeWithVertex[1][1]=n4;
        edgeWithVertex[2][0]=n2; edgeWithVertex[2][1]=n4;
        edgeWithVertex[3][0]=n2; edgeWithVertex[3][1]=n3;
    }
};

//------------------------------------------------------------------------------

int getPolygons(Elem &elem, LevelSetStructure &LSS, PolygonReconstructionData* polygons);
int getPolygon(Elem &elem, LevelSetStructure &LSS, int polygon[4][2]);
void getPolygonNormal(SVec<double,3>& X, Vec3D &normal, LevelSetStructure &LSS, PolygonReconstructionData &polygon);

#endif
