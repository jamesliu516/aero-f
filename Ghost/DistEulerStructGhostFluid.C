#include <Domain.h>
#include "Ghost/DistEulerStructGhostFluid.h"

//-----------------------------------------------------------------------------------

template<int dim>
void DistEulerStructGhostFluid::setupCommunication(Domain *domain, DistSVec<double,3>*X, DistSVec<double,dim> &U)
{  // At this moment doesn't need U, but might need it in future.
    if (solidsurface)  prepareForCommunication();
    else {com->fprintf(stderr,"ERROR: Solid surface file not specified. Aborting...\n"); exit(-1);}
    if (!givenBB) specifyBoundingBox(X);
    int m, n, mn;
//    double dx = specifydx(domain,X); 
//    m = static_cast<int>((Xmax - Xmin) / dx);
//    n = static_cast<int>((Ymax - Ymin) / dx);
//    mn = static_cast<int>((Zmax - Zmin) / dx);
    m = 256; n = mn = 86; //for debug only. should be specified by "specifydx(...)"
    com->fprintf(stderr,"For Cart. Mesh, dx = %lf, dy = %lf, dz = %lf.\n", (Xmax-Xmin)/m, (Ymax-Ymin)/n, (Zmax-Zmin)/mn);
    //specifyBandwidth();  
    //getTetNearInterface(*X, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax);
    initializePhysBAM(m,n,mn,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax);
    //initializePhysBAMMPI();

    computeLevelSet();  //compute level-set on Cart. grid.

    getPhiFromModule(*X,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax); //interpolate Phi onto fluid (tet) grid.
}

//----------------------------------------------------------------------------------------

