//#####################################################################
// Copyright 2007-2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Continuous_Collision_Detection/POINT_FACE_COLLISION.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include "PHYSBAM_INTERFACE.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> PhysBAMInterface<T>::
PhysBAMInterface(TRIANGLE_MESH& triangle_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input)
    :thickness_parameter(1e-6),thickness_over_two(5e-7),triangle_mesh(triangle_mesh_input),
    triangle_list(triangle_mesh.elements.m),triangle_hierarchy(0),SubD(0),particles(particles_input)
{
    Update(0,true);

    saved_state.y=new GEOMETRY_PARTICLES<TV>;
    saved_state.y->array_collection.Add_Elements(particles.array_collection.Size());
    saved_state.x=(T)0;saved_state.y->X=particles.X;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> PhysBAMInterface<T>::
~PhysBAMInterface()
{
    delete &triangle_mesh;
    delete triangle_hierarchy;
    delete &particles;
    delete saved_state.y;
    for(int iSub=1;iSub<=SubD.Size();++iSub){
        delete SubD(iSub).scoped_triangle_mesh;
        delete SubD(iSub).triangle_hierarchy;}
}
//#####################################################################
// Update
//#####################################################################
template<class T> void PhysBAMInterface<T>::
Update(const int numLocSub,const bool rebuild_hierarchy)
{ // Only updates the GLOBAL mesh and hierarchy.
    if(SubD.Size() < numLocSub) SubD.Resize(numLocSub);
    for(int t=1;t<=triangle_mesh.elements.m;t++)
        triangle_list(t)=TRIANGLE_3D<T>(particles.X.Subset(triangle_mesh.elements(t)));

    if(rebuild_hierarchy){
        delete triangle_hierarchy;triangle_hierarchy=new TRIANGLE_HIERARCHY<T>(triangle_mesh,particles,triangle_list,true,0);}
    else triangle_hierarchy->Update_Boxes();
}
//#####################################################################
// UpdateScope
//#####################################################################
template<class T> void PhysBAMInterface<T>::
UpdateScope(const int subD,const bool use_global_scope)
{ // Updates the SCOPED mesh and hierarchy, and cleans out acceleration structures.
    SubDInterface<T>& sub = SubD(subD);

    if(use_global_scope){
        sub.scoped_triangle_mesh=&triangle_mesh;
        sub.triangle_list=triangle_list;
        sub.triangle_hierarchy=triangle_hierarchy;}
    else{
        sub.scope.Resize(0);
        for(std::set<int>::const_iterator iter=sub.next_scope.begin();iter!=sub.next_scope.end();iter++) sub.scope.Append(*iter);

        ARRAY<VECTOR<int,3> > scoped_elements(sub.scope.m);
        for(int t=1;t<=scoped_elements.m;++t) scoped_elements(t)=triangle_mesh.elements(sub.scope(t));
        delete sub.scoped_triangle_mesh; sub.scoped_triangle_mesh=new TRIANGLE_MESH(triangle_mesh.number_nodes,scoped_elements);

        sub.triangle_list.Resize(sub.scope.m);
        for(int t=1;t<=sub.scope.m;t++)
            sub.triangle_list(t)=TRIANGLE_3D<T>(particles.X.Subset(sub.scoped_triangle_mesh->elements(t)));
        delete sub.triangle_hierarchy;sub.triangle_hierarchy=new TRIANGLE_HIERARCHY<T>(*sub.scoped_triangle_mesh,particles,sub.triangle_list,true,0);}


    sub.candidates.Resize(0);
    sub.next_scope.clear();
}
//#####################################################################
// SaveOldState
//#####################################################################
template<class T> void PhysBAMInterface<T>::
SaveOldState()
{
    saved_state.x+=(T)1;saved_state.y->X=particles.X;
}
//#####################################################################
// HasCloseTriangle
//#####################################################################
template<class T> bool PhysBAMInterface<T>::
HasCloseTriangle(const int subD,const TV position,const TV min_corner,const TV max_corner,int* index,bool* is_occluded,ARRAY<int>* cand)
{
    SubDInterface<T>& sub = SubD(subD);
    ARRAY<int> candidates;
    RANGE<TV> bounding_box(min_corner, max_corner);
    sub.triangle_hierarchy->Intersection_List(bounding_box,candidates,thickness_over_two);

    if(candidates.Size()>0){
        if(index) *index = sub.candidates.Append(candidates);
        for(int i=1;i<=candidates.Size();++i) sub.next_scope.insert(sub.scope(candidates(i)));
        if(cand) for(int i=1;i<=candidates.Size();++i) cand->Append(sub.scope(candidates(i)));
        if(is_occluded){*is_occluded=false;
            for(int t=1;t<=candidates.Size() && !*is_occluded;++t) *is_occluded=(sub.triangle_list(candidates(t)).Point_Inside_Triangle(position,thickness_over_two));}
        return true;}
    else return false;
}
//#####################################################################
// Intersect
//#####################################################################
template<class T> IntersectionResult<T> PhysBAMInterface<T>::
Intersect(const TV& start,const TV& end,const T thickness) const
{
    IntersectionResult<T> result;
    RAY<TV> edge(SEGMENT_3D<T>(start,end));
    TV weights;T old_t_max=edge.t_max;
    if(triangle_hierarchy->Intersection(edge,weights,thickness*(T).5,true)){ // TODO: Pull out to avoid hierarchy traversal
        result.triangleID=edge.aggregate_id;
        result.zeta[0]=weights(1);result.zeta[1]=weights(2);result.zeta[2]=weights(3);
        result.alpha=(T)1-edge.t_max/old_t_max;}
    else
        result.triangleID=-1;
    return result;
}
//#####################################################################
// Intersect
//#####################################################################
namespace{
template<class T> void
Assign_Intersection_Information(const ARRAY<int>& scope,const RAY<VECTOR<T,3> >& intersection_ray,const RAY<VECTOR<T,3> >& original_ray,const VECTOR<T,3>& weights,
                                IntersectionResult<T>& result)
{
    if(intersection_ray.intersection_location == RAY<VECTOR<T,3> >::INTERIOR_POINT){
        result.triangleID=scope(intersection_ray.aggregate_id);
        result.alpha=(T)1-intersection_ray.t_max/original_ray.t_max;
        result.zeta[0]=weights(1);result.zeta[1]=weights(2);result.zeta[2]=weights(3);}
    else result.triangleID=-1;
}

template<class T> void
Assign_Intersection_Information(const int triangleID,const T alpha,const VECTOR<T,3>& weights,
                                IntersectionResult<T>& result)
{
    result.triangleID=triangleID;
    result.alpha=alpha;
    result.zeta[0]=weights(1);result.zeta[1]=weights(2);result.zeta[2]=weights(3);
}

template<class T> void
Ray_Intersection(const ARRAY<TRIANGLE_3D<T> >& triangle_list,const ARRAY<int>& candidates,const T thickness_over_two,RAY<VECTOR<T,3> >& ray,VECTOR<T,3>& weights)
{
    bool intersected=false;
    for(int t=1;t<=candidates.m;++t){
        RAY<VECTOR<T,3> > ray_temp=ray;
        if(INTERSECTION::Intersects(ray,triangle_list(candidates(t)),thickness_over_two)) {
            ray.aggregate_id=candidates(t);intersected=true;}}
    if(intersected) triangle_list(ray.aggregate_id).Closest_Point(ray.Point(ray.t_max),weights);
}

template<class T> bool
Point_Inside_Triangle(const ARRAY<TRIANGLE_3D<T> >& triangle_list,const VECTOR<T,3>& point,const T thickness_over_two, int& triangleID, VECTOR<T,3>& weights,
                      const ARRAY<int>* candidates=0)
{
    if(candidates){
        for(int i=1;i<=candidates->Size();++i) if(triangle_list((*candidates)(i)).Point_Inside_Triangle(point,thickness_over_two)){
            triangleID=(*candidates)(i);weights=triangle_list((*candidates)(i)).Barycentric_Coordinates(point);return true;}}
    else {
        for(int i=1;i<=triangle_list.Size();++i) if(triangle_list(i).Point_Inside_Triangle(point,thickness_over_two)){
            triangleID=i;weights=triangle_list(i).Barycentric_Coordinates(point);return true;}}
    return false;
}

template<class T> int
Intersection_Helper(const ARRAY<TRIANGLE_3D<T> >& triangle_list,const ARRAY<int>& scope,const ARRAY<int>& candidates,
                    const VECTOR<T,3>& start,const VECTOR<T,3>& end,const bool start_node_occluded,const bool end_node_occluded,
                    const T thickness_over_two,IntersectionResult<T>& result){
    const int PERMITTED_ATTEMPTS=3;
    int retryAttempts=0;
    T modified_thickness_over_two=thickness_over_two;
    typedef VECTOR<T,3> TV;
    RAY<TV> intersection_ray(SEGMENT_3D<T>(start,end)),original_ray(SEGMENT_3D<T>(start,end));
    TV weights(0,0,0);

    if(!start_node_occluded) // The usual case
        Ray_Intersection(triangle_list,candidates,thickness_over_two,intersection_ray,weights);
    else{ // if(start_node_occluded)
        // LOG::cerr<<"("<<start.x<<", "<<start.y<<", "<<start.z<<")"<<" -> "<<"("<<end.x<<", "<<end.y<<", "<<end.z<<")"<<": STARTING (start_node_occluded)"<<std::endl;
        int triangleID=-10;
        modified_thickness_over_two=thickness_over_two;
        for(int number_of_attempts=0;number_of_attempts<PERMITTED_ATTEMPTS;++number_of_attempts) {
            if(number_of_attempts) ++retryAttempts;
            if(Point_Inside_Triangle(triangle_list,start,modified_thickness_over_two,triangleID,weights,&candidates)) {
                // LOG::cerr<<"Point_Inside_Triangle 1 returned triangleID = "<<triangleID<<std::endl;
                Assign_Intersection_Information(scope(triangleID),modified_thickness_over_two,weights,result);break;}
            modified_thickness_over_two*=(T)2;}
        // LOG::cerr<<"("<<start.x<<", "<<start.y<<", "<<start.z<<")"<<" -> "<<"("<<end.x<<", "<<end.y<<", "<<end.z<<")"<<": (start_node_occluded) PART 1, "<<triangleID<<std::endl;
        if(triangleID <= 0 && Point_Inside_Triangle(triangle_list,start,modified_thickness_over_two,triangleID,weights)) {
                // LOG::cerr<<"Point_Inside_Triangle 2 returned triangleID = "<<triangleID<<std::endl;
            Assign_Intersection_Information(scope(triangleID),modified_thickness_over_two,weights,result);}
        // LOG::cerr<<"("<<start.x<<", "<<start.y<<", "<<start.z<<")"<<" -> "<<"("<<end.x<<", "<<end.y<<", "<<end.z<<")"<<": (start_node_occluded) PART 2, "<<triangleID<<std::endl;
        if(triangleID <= 0){
            LOG::cerr<<"An occluded node is found that appears to be visible!"<<std::endl;Assign_Intersection_Information(-1,(T)0,weights,result);}
        return retryAttempts;}

    if(intersection_ray.intersection_location == RAY<TV>::START_POINT) {
        // LOG::cerr<<"EDGE case A occuring at "<<"("<<start.x<<", "<<start.y<<", "<<start.z<<")"<<" -> "<<"("<<end.x<<", "<<end.y<<", "<<end.z<<")"<<std::endl;
        for(int number_of_attempts=0;number_of_attempts<PERMITTED_ATTEMPTS && intersection_ray.intersection_location == RAY<TV>::START_POINT;++number_of_attempts) {
            ++retryAttempts;
            modified_thickness_over_two=thickness_over_two*(T).5;
            intersection_ray.Restore_Intersection_Information(original_ray);
            Ray_Intersection(triangle_list,candidates,modified_thickness_over_two,intersection_ray,weights);}
        if(intersection_ray.intersection_location == RAY<TV>::START_POINT){LOG::cerr<<"A visible node appears to be occluded!"<<std::endl;Assign_Intersection_Information(-2,(T)0,weights,result);return retryAttempts;}}

    if(end_node_occluded && intersection_ray.intersection_location != RAY<TV>::INTERIOR_POINT) {
        // LOG::cerr<<"EDGE case B occuring at "<<"("<<start.x<<", "<<start.y<<", "<<start.z<<")"<<" -> "<<"("<<end.x<<", "<<end.y<<", "<<end.z<<")"<<std::endl;
        modified_thickness_over_two=thickness_over_two;
        int triangleID=-2;
        for(int number_of_attempts=0;number_of_attempts<PERMITTED_ATTEMPTS;++number_of_attempts) {
            ++retryAttempts;
            if(Point_Inside_Triangle(triangle_list,end,modified_thickness_over_two,triangleID,weights,&candidates)) {
                // LOG::cerr<<"Point_Inside_Triangle 3 returned triangleID = "<<triangleID<<std::endl;
                Assign_Intersection_Information(scope(triangleID),(T)1-modified_thickness_over_two,weights,result);break;}
            modified_thickness_over_two*=(T)2;}
        if(triangleID <= 0 && Point_Inside_Triangle(triangle_list,end,modified_thickness_over_two,triangleID,weights)) {
                // LOG::cerr<<"Point_Inside_Triangle 4 returned triangleID = "<<triangleID<<std::endl;
            Assign_Intersection_Information(scope(triangleID),(T)1-modified_thickness_over_two,weights,result);}
        if(triangleID <= 0){LOG::cerr<<"An occluded node is reported visible!"<<std::endl;Assign_Intersection_Information(-3,(T)0,weights,result);}}
    else Assign_Intersection_Information(scope,intersection_ray,original_ray,weights,result);

    return retryAttempts;
}
}

template<class T> void PhysBAMInterface<T>::
Intersect(const int subD,const ARRAY<TV>& node_positions,const ARRAY<bool>& occluded_node,ARRAY<TRIPLE<VECTOR<int,3>,IntersectionResult<T>,IntersectionResult<T> > >& edges_and_results) const
{
    const SubDInterface<T>& sub = SubD(subD);

    int retryAttempts=0;
    for(int i=1;i<=edges_and_results.m;++i){
        VECTOR<int,3>& edge_nodes=edges_and_results(i).x;
        int start_node=edge_nodes(1),end_node=edge_nodes(2);
        VECTOR<T,3> left_node=node_positions(start_node),right_node=node_positions(end_node);
        bool left_node_occluded=occluded_node(start_node),right_node_occluded=occluded_node(end_node);
        ARRAY<int> candidates; // Perform an intersect on candidate lists.
        for(int j=1;j<=sub.candidates(start_node).Size();++j) for(int k=1;k<=sub.candidates(end_node).Size();++k)
            if(sub.candidates(start_node)(j)==sub.candidates(end_node)(k)) candidates.Append(sub.candidates(start_node)(j));

        retryAttempts += Intersection_Helper(sub.triangle_list,sub.scope,candidates,left_node,right_node,
                                             left_node_occluded,right_node_occluded,thickness_over_two,edges_and_results(i).y);
        retryAttempts += Intersection_Helper(sub.triangle_list,sub.scope,candidates,right_node,left_node,
                                             right_node_occluded,left_node_occluded,thickness_over_two,edges_and_results(i).z);}
}
//#####################################################################
// computeSweptNodes
//#####################################################################
template<class T> void PhysBAMInterface<T>::
computeSweptNodes(const int subD,const ARRAY<TV>& node_positions,ARRAY<bool>& swept_node,const T dt) const
{ // TODO(jontg): Assuming a 'normalized' dt of 1 for now... verify that this is an OK assumption to make!
    const SubDInterface<T>& sub = SubD(subD);
    T relative_speed,collision_time;VECTOR<T,3> normal,weights; // Unused
    ARRAYS_COMPUTATIONS::Fill(swept_node,false);
    for(int i=1;i<=node_positions.Size();++i){
        const ARRAY<int>& candidates(sub.candidates(i));
        for(int t=1;t<=candidates.Size() && !swept_node(i);++t){
            const TRIANGLE_3D<T> initial_simplex(saved_state.y->X.Subset(sub.scoped_triangle_mesh->elements(candidates(t))));
            const TRIANGLE_3D<T>& terminal_simplex(sub.triangle_list(candidates(t)));
            swept_node(i) = (POINT_SIMPLEX_NO_COLLISION != CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Robust_Point_Triangle_Collision(initial_simplex,terminal_simplex,
                                                           node_positions(i),node_positions(i),dt,thickness_over_two,collision_time,normal,weights,relative_speed));}}
}
//#####################################################################
template class PhysBAMInterface<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PhysBAMInterface<double>;
#endif
