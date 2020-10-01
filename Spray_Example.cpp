//!#####################################################################
//! \file Spray_Example.cpp
//!#####################################################################
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Initializer.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Iterator.h>
#include <nova/Dynamics/Hierarchy/Advection/Grid_Hierarchy_Advection.h>
#include <nova/Tools/Grids/Grid_Iterator_Face.h>
#include <nova/Tools/Krylov_Solvers/Conjugate_Gradient.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include "Boundary_Value_Initializer.h"
#include "Uniform_Velocity_Field_Initializer.h"
#include "Velocity_Field_Traverser.h"
#include "Boundary_Condition_Helper.h"
#include "Compute_Time_Step.h"
#include "Density_Modifier.h"
#include "Source_Velocity_Setup.h"
#include "Source_Adder.h"
#include "Initialize_Dirichlet_Cells.h"
#include "Spray_Example.h"
#include "Ficks_RHS_Helper.h"
#include "Lap_Calculator.h"
#include "Density_Backup_Helper.h"
#include "Explicit_Face_Qc_Updater.h"
#include "Implicit_Face_Qc_Updater.h"
#include "Uniform_Grid_Helper/Uniform_Grid_Advection_Helper.h"
#include "Uniform_Grid_Helper/Uniform_Grid_Averaging_Helper.h"
#include "Density_Clamp_Helper.h"
#include "Flip_Helper.h"
#include "Write_To_File_Helper.h"
#include <omp.h>
#include <chrono>

#include "Alpha_Helper.h"
#include "Clamp_Alpha_Helper.h"
#include "Relative_Velocity_Helper.h"
#include "Compute_Drag_Force_Helper.h"
#include "Apply_Drag_Force_Helper.h"
#include "Combination_RHS_Helper.h"
#include "Poisson_Solver/Poisson_CG_System.h"
#include "Apply_Combination_Pressure.h"
#include "Apply_External_Force_Helper.h"
#include "Check_Velocity_Helper.h"
#include "Exp_Helper.h"

using namespace std::chrono;
using namespace Nova;
namespace Nova{
extern int number_of_threads;
}
//######################################################################
// Constructor
//######################################################################
template<class T,int d> Spray_Example<T,d>::
Spray_Example()
    :Base(),hierarchy(nullptr),rasterizer(nullptr)
{
    level=0;
    init_v1=TV({0.,-1.});   init_v2=TV({-3.,0.});
    face_velocity1_channels(0)                          = &Struct_type::ch0;
    face_velocity1_channels(1)                          = &Struct_type::ch1;
    if(d==3) face_velocity1_channels(2)                 = &Struct_type::ch2;
    face_velocity1_backup_channels(0)                   = &Struct_type::ch3;
    face_velocity1_backup_channels(1)                   = &Struct_type::ch4;
    if(d==3) face_velocity1_backup_channels(2)          = &Struct_type::ch5;

    face_velocity2_channels(0)                          = &Struct_type::ch6;
    face_velocity2_channels(1)                          = &Struct_type::ch7;
    if(d==3) face_velocity2_channels(2)                 = &Struct_type::ch8;
    face_velocity2_backup_channels(0)                   = &Struct_type::ch9;
    face_velocity2_backup_channels(1)                   = &Struct_type::ch10;
    if(d==3) face_velocity2_backup_channels(2)          = &Struct_type::ch11;

    alpha1_channel                                      = &Struct_type::ch12;
    alpha2_channel                                      = &Struct_type::ch13;
}
//######################################################################
// Initialize
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Initialize()
{
    diffusion_rt=(T)0.; qc_advection_rt=(T)0.; 
    Initialize_SPGrid();
    Initialize_Fluid_State(test_number);

    const Grid<T,d>& grid=hierarchy->Lattice(level);
    Log::cout<<"dx: "<<grid.dX<<std::endl;
    Log::cout<<"resolution: "<<grid.counts<<std::endl;
    Log::cout<<"domain: "<<grid.domain.min_corner<<", "<<grid.domain.max_corner<<std::endl;
}
//######################################################################
// Initialize_SPGrid
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Initialize_SPGrid()
{
    Log::Scope scope("Initialize_SPGrid");
    Initialize_Rasterizer(test_number);
    for(Grid_Hierarchy_Iterator<d,Hierarchy_Rasterizer> iterator(hierarchy->Lattice(levels-1).Cell_Indices(),levels-1,*rasterizer);iterator.Valid();iterator.Next());
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Ghost_Cells(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Valid_Faces(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Active_Faces(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Active_Nodes(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Shared_Nodes(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_Ghost_Nodes(*hierarchy);
    Grid_Hierarchy_Initializer<Struct_type,T,d>::Flag_T_Junction_Nodes(*hierarchy);
    Initialize_Dirichlet_Cells<Struct_type,T,d>(*hierarchy,domain_walls);
    Set_Neumann_Faces_Inside_Sources();
    hierarchy->Update_Block_Offsets();
    hierarchy->Initialize_Red_Black_Partition(2*number_of_threads);
}
//######################################################################
// Limit_Dt
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Limit_Dt(T& dt,const T time)
{
    T dt_convection=(T)0.;

    Vector<uint64_t,d> other_face_offsets;
    for(int axis=0;axis<d;++axis) other_face_offsets(axis)=Topology_Helper::Axis_Vector_Offset(axis);

    // doubt
    Compute_Time_Step<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity1_channels,
                                        other_face_offsets,level,dt_convection);

    if(dt_convection>(T)1e-5) dt=cfl/dt_convection;
    Log::cout<<"Time Step: "<<dt<<std::endl;
}
//######################################################################
// Advect_Alpha
//######################################################################
template    <class T,int d> void Spray_Example<T,d>::
Advect_Alpha(const T& dt)
{
    Log::cout<<"advecting alpha ..."<<std::endl;
    Channel_Vector cell_velocity_channels;
    cell_velocity_channels(0)               = &Struct_type::ch14;
    cell_velocity_channels(1)               = &Struct_type::ch15;
    if(d==3) cell_velocity_channels(2)      = &Struct_type::ch16;
    T Struct_type::* temp_channel           = &Struct_type::ch17;
    Vector<uint64_t,d> other_face_offsets;
    for(int axis=0;axis<d;++axis) other_face_offsets(axis)=Topology_Helper::Axis_Vector_Offset(axis);
    Uniform_Grid_Averaging_Helper<Struct_type,T,d>::Uniform_Grid_Average_Face_Velocities_To_Cells(*hierarchy,hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity1_channels,cell_velocity_channels,other_face_offsets);
    Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Density(*hierarchy,cell_velocity_channels,alpha1_channel,temp_channel,dt);
    // ensure alpha_air + alpha_water = 1
    Alpha_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),alpha1_channel,alpha2_channel);
    Clamp_Alpha_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),alpha1_channel,alpha2_channel);
}
//######################################################################
// Backup_Velocity
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Backup_Velocity()
{
    Log::cout<<"backing up velocity ..."<<std::endl;
    for(int axis=0;axis<d;++axis){ 
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity1_backup_channels(axis));
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity2_backup_channels(axis));
        SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity1_channels(axis),
                                        face_velocity1_backup_channels(axis),(unsigned)Topology_Helper::Face_Active_Mask(axis));
        SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity2_channels(axis),
                                        face_velocity2_backup_channels(axis),(unsigned)Topology_Helper::Face_Active_Mask(axis));}
    // Check_Velocity_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity1_channels,face_velocity2_channels);
}
//######################################################################
// Ficks_Diffusion
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Apply_Drag_Force(const T& dt)
{
    Log::cout<<"applying drag force ..."<<std::endl;
    // enum {number_of_nodes_per_cell                              = Topology_Helper::number_of_nodes_per_cell};
    // enum {number_of_nodes_per_face                              = Topology_Helper::number_of_nodes_per_face};
    // Channel_Vector interpolated_rel_face_velocity_channels;
    // interpolated_rel_face_velocity_channels(0)                  = &Struct_type::ch14;
    // interpolated_rel_face_velocity_channels(1)                  = &Struct_type::ch15;
    // if(d==3) interpolated_rel_face_velocity_channels(2)         = &Struct_type::ch16;
    // Channel_Vector rel_face_velocity_channels;
    // rel_face_velocity_channels(0)                               = &Struct_type::ch17;
    // rel_face_velocity_channels(1)                               = &Struct_type::ch18;
    // if(d==3) rel_face_velocity_channels(2)                      = &Struct_type::ch19;
    // Channel_Vector drag_force_channels;
    // drag_force_channels(0)                                      = &Struct_type::ch20;
    // drag_force_channels(1)                                      = &Struct_type::ch21;
    // if(d==3) drag_force_channels(2)                             = &Struct_type::ch22;
    // Vector<uint64_t,d> negative_face_offsets; for(int axis=0;axis<d;++axis) negative_face_offsets(axis)=Topology_Helper::Negative_Axis_Vector_Offset(axis);
    // uint64_t nodes_of_cell_offsets[number_of_nodes_per_cell]; Topology_Helper::Nodes_Of_Cell_Offsets(nodes_of_cell_offsets);
    // for(int axis=0;axis<d;++axis) { 
    //     SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),rel_face_velocity_channels(axis));
    //     SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),drag_force_channels(axis));
    //     // false: 
    //     Relative_Velocity_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity1_channels,face_velocity2_channels,rel_face_velocity_channels,true);
    //     // interpolate velocity at face
    //     Uniform_Grid_Averaging_Helper<Struct_type,T,d>::Uniform_Grid_Average_Face_Velocities_To_Faces(*hierarchy,hierarchy->Allocator(level),hierarchy->Blocks(level),
    //                                                 rel_face_velocity_channels,interpolated_rel_face_velocity_channels,negative_face_offsets,nodes_of_cell_offsets,axis);
    //     Compute_Drag_Force_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),drag_force_channels(axis),alpha1_channel,interpolated_rel_face_velocity_channels,rho2,axis);
    // }  

    for(int axis=0;axis<d;++axis)
        Exp_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity1_channels(axis),face_velocity2_channels(axis),alpha1_channel,
                                    dt,rho1,rho2,axis);
        // Apply_Drag_Force_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity1_channels(axis),face_velocity2_channels(axis),drag_force_channels(axis),alpha1_channel,alpha2_channel,dt,rho1,rho2,axis);
    
}
//######################################################################
// Apply_External_Force
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Apply_External_Force(const T& dt)
{
    Log::cout<<"applying external force ..."<<std::endl;
    TV gravity=TV::Axis_Vector(1)*(T)-9.8;
    Apply_External_Force_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity1_channels,gravity,dt);
    Apply_External_Force_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity2_channels,gravity,dt);
}
//######################################################################
// Combination_Project
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Combination_Project()
{
    Log::cout<<"projecting ..."<<std::endl;
    // set up divergence channel
    T Struct_type::* rhs_channel            = &Struct_type::ch14;
    T Struct_type::* pressure_channel       = &Struct_type::ch15;

    Log::cout<<"before clearing ..."<<std::endl;
    // Check_Velocity_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity1_channels,face_velocity2_channels);


    // // clear
    // for(int level=0;level<levels;++level){
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),pressure_channel);
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),rhs_channel);

    const T one_over_dx=hierarchy->Lattice(level).one_over_dX(0);
    Array<TV> source_velocity1;         Array<TV> source_velocity2;
    source_velocity1.Append(init_v1);   source_velocity2.Append(init_v2);
    
    Log::cout<<"before setting boundary ..."<<std::endl;
    // Check_Velocity_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity1_channels,face_velocity2_channels);

    // set boundary conditions
    Boundary_Condition_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity1_channels,pressure_channel,level);
    Boundary_Condition_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity2_channels,pressure_channel,level);
    Source_Velocity_Setup<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),velocity1_sources,source_velocity1,face_velocity1_channels,level);
    Source_Velocity_Setup<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),velocity2_sources,source_velocity2,face_velocity2_channels,level);

    Log::cout<<"before projecting ..."<<std::endl;
    // Check_Velocity_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity1_channels,face_velocity2_channels);


    Combination_RHS_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Allocator(level),hierarchy->Blocks(level),rhs_channel,alpha1_channel,alpha2_channel,face_velocity1_channels,face_velocity2_channels,rho1,rho2,one_over_dx);

    Poisson_CG_System<Struct_type,T,d> cg_system(*hierarchy,alpha1_channel,rho1,rho2);

    T Struct_type::* q_channel              = &Struct_type::ch16;
    T Struct_type::* r_channel              = &Struct_type::ch17;
    T Struct_type::* s_channel              = &Struct_type::ch18;
    T Struct_type::* k_channel              = &Struct_type::ch19;
    T Struct_type::* z_channel              = &Struct_type::ch19;

    // clear all channels
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),q_channel);
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),r_channel);
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),s_channel);
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),z_channel);

    CG_Vector<Struct_type,T,d> x_V(*hierarchy,pressure_channel);
    CG_Vector<Struct_type,T,d> b_V(*hierarchy,rhs_channel);
    CG_Vector<Struct_type,T,d> q_V(*hierarchy,q_channel);
    CG_Vector<Struct_type,T,d> r_V(*hierarchy,r_channel);
    CG_Vector<Struct_type,T,d> s_V(*hierarchy,s_channel);
    CG_Vector<Struct_type,T,d> k_V(*hierarchy,k_channel);
    CG_Vector<Struct_type,T,d> z_V(*hierarchy,z_channel);

    Conjugate_Gradient<T> cg;
    cg.print_diagnostics=false; cg.print_residuals=false;
    cg_system.Multiply(x_V,r_V);
    r_V-=b_V;
    const T b_norm=cg_system.Convergence_Norm(r_V);
    Log::cout<<"Norm: "<<b_norm<<std::endl;
    cg.restart_iterations=cg_restart_iterations;
    const T tolerance=std::max((T)1e-10*b_norm,(T)1e-10);
    cg.Solve(cg_system,x_V,b_V,q_V,s_V,r_V,k_V,z_V,tolerance,0,cg_iterations);

    Apply_Combination_Pressure<Struct_type,T,d>(*hierarchy,hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity2_channels,pressure_channel,rho1,one_over_dx);
    Apply_Combination_Pressure<Struct_type,T,d>(*hierarchy,hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity1_channels,pressure_channel,rho2,one_over_dx);
    // Check_Velocity_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity1_channels,face_velocity2_channels);
}
//######################################################################
// Add_Source
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Add_Source(const T& dt)
{
    Log::cout<<"adding source ..."<<std::endl;
    Source_Adder<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),alpha1_channel,density_sources,source_rate,dt,level);
    Alpha_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),alpha1_channel,alpha2_channel);
    Clamp_Alpha_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),alpha1_channel,alpha2_channel);
}
//######################################################################
// Advect_Face_Velocities
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Advect_Face_Velocities(const T& dt)
{
    Log::cout<<"advecting face velocities ..."<<std::endl;
    Channel_Vector interpolated_face_velocity_channels;
    interpolated_face_velocity_channels(0)              = &Struct_type::ch14;
    interpolated_face_velocity_channels(1)              = &Struct_type::ch15;
    if(d==3) interpolated_face_velocity_channels(2)     = &Struct_type::ch16;
    Channel_Vector face_velocity_backup_channels;
    face_velocity_backup_channels(0)                    = &Struct_type::ch17;
    face_velocity_backup_channels(1)                    = &Struct_type::ch18;
    if(d==3) face_velocity_backup_channels(2)           = &Struct_type::ch19;
    T Struct_type::* temp_channel                       = &Struct_type::ch20;
    // phase 1
    // backup face velocity
    for(int axis=0;axis<d;++axis) { SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity_backup_channels(axis));
        SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity1_channels(axis),
                                            face_velocity_backup_channels(axis),Topology_Helper::Face_Active_Mask(axis));}  
    Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Face_Velocities(*hierarchy,face_velocity1_channels,face_velocity_backup_channels,interpolated_face_velocity_channels,temp_channel,dt);
    // phase 2
    // backup face velocity
    for(int axis=0;axis<d;++axis) { SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity_backup_channels(axis));
        SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity2_channels(axis),
                                            face_velocity_backup_channels(axis),Topology_Helper::Face_Active_Mask(axis));}  
    Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Face_Velocities(*hierarchy,face_velocity2_channels,face_velocity_backup_channels,interpolated_face_velocity_channels,temp_channel,dt);

    // Check_Velocity_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity1_channels,face_velocity2_channels);
}
//######################################################################
// Set_Neumann_Faces_Inside_Sources
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Set_Neumann_Faces_Inside_Sources()
{
    auto flags=hierarchy->Allocator(level).template Get_Array<Struct_type,unsigned>(&Struct_type::flags);
    for(Grid_Iterator_Face<T,d> iterator(hierarchy->Lattice(level));iterator.Valid();iterator.Next()){
        const int axis=iterator.Axis();const T_INDEX& face_index=iterator.Face_Index();
        uint64_t face_offset=Flag_array_mask::Linear_Offset(face_index._data);
        const uint64_t face_active_mask=Topology_Helper::Face_Active_Mask(axis);
        const uint64_t face_valid_mask=Topology_Helper::Face_Valid_Mask(axis);
        const TV X=hierarchy->Lattice(level).Face(axis,face_index);
        for(size_t i=0;i<velocity1_sources.size();++i) if(velocity1_sources(i)->Inside(X)||velocity2_sources(i)->Inside(X)){
            if(hierarchy->template Set<unsigned>(level,&Struct_type::flags).Is_Set(face_offset,face_valid_mask))
                flags(face_offset)&=~face_active_mask;break;}}
}
//######################################################################
// Initialize_Velocity_Field
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Initialize_Velocity_Field()
{
    Log::cout<<"initializing velocity_field ..."<<std::endl;
    Array<TV> source_velocity1; Array<TV> source_velocity2;
    source_velocity1.Append(init_v1); source_velocity2.Append(init_v2);
    Source_Velocity_Setup<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),velocity1_sources,source_velocity1,face_velocity1_channels,level);
    Source_Velocity_Setup<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),velocity2_sources,source_velocity2,face_velocity2_channels,level);
    // Check_Velocity_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity1_channels,face_velocity2_channels);
}
//######################################################################
// Register_Options
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Register_Options()
{
    Base::Register_Options();

    parse_args->Add_Double_Argument("-cfl",(T).5,"CFL number.");
    parse_args->Add_Integer_Argument("-levels",1,"Number of levels in the SPGrid hierarchy.");
    parse_args->Add_Integer_Argument("-test_number",1,"Test number.");
    parse_args->Add_Integer_Argument("-mg_levels",1,"Number of levels in the Multigrid hierarchy.");
    parse_args->Add_Integer_Argument("-threads",1,"Number of threads for OpenMP to use");
    if(d==2) parse_args->Add_Vector_2D_Argument("-size",Vector<double,2>(64.),"n","Grid resolution");
    else if(d==3) parse_args->Add_Vector_3D_Argument("-size",Vector<double,3>(64.),"n","Grid resolution");
    if(d==2) parse_args->Add_Vector_2D_Argument("-domain",Vector<double,2>(1.),"n","domain");
    else if(d==3) parse_args->Add_Vector_3D_Argument("-domain",Vector<double,3>(1.),"n","domain");
    parse_args->Add_Double_Argument("-rho1",(T)10.,"rho1.");
    parse_args->Add_Double_Argument("-rho2",(T)1.,"rho2.");
    parse_args->Add_Double_Argument("-sr",(T)1.,"Source rate");
    parse_args->Add_Option_Argument("-uvf","Uniform velocity field");
    // for CG
    parse_args->Add_Integer_Argument("-cg_iterations",100,"Number of CG iterations.");
    parse_args->Add_Integer_Argument("-cg_restart_iterations",40,"Number of CG restart iterations.");
    parse_args->Add_Double_Argument("-cg_tolerance",1e-4,"CG tolerance");
}
//######################################################################
// Parse_Options
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Parse_Options()
{
    Base::Parse_Options();

    cfl=(T)parse_args->Get_Double_Value("-cfl");
    levels=parse_args->Get_Integer_Value("-levels");
    test_number=parse_args->Get_Integer_Value("-test_number");
    mg_levels=parse_args->Get_Integer_Value("-mg_levels");
    number_of_threads=parse_args->Get_Integer_Value("-threads");
    omp_set_num_threads(number_of_threads);
    if(d==2){auto cell_counts_2d=parse_args->Get_Vector_2D_Value("-size");for(int v=0;v<d;++v) counts(v)=cell_counts_2d(v);}
    else{auto cell_counts_3d=parse_args->Get_Vector_3D_Value("-size");for(int v=0;v<d;++v) counts(v)=cell_counts_3d(v);}
    if(d==2){auto cell_domain_2d=parse_args->Get_Vector_2D_Value("-domain");for(int v=0;v<d;++v) domain(v)=cell_domain_2d(v);}
    else{auto cell_domain_3d=parse_args->Get_Vector_3D_Value("-domain");for(int v=0;v<d;++v) domain(v)=cell_domain_3d(v);}
    rho1=parse_args->Get_Double_Value("-rho1");
    rho2=parse_args->Get_Double_Value("-rho2");
    source_rate=parse_args->Get_Double_Value("-sr");
    uvf=parse_args->Get_Option_Value("-uvf");
    cg_iterations=parse_args->Get_Integer_Value("-cg_iterations");
    cg_restart_iterations=parse_args->Get_Integer_Value("-cg_restart_iterations");
    cg_tolerance=(T)parse_args->Get_Double_Value("-cg_tolerance");
    switch (test_number){
    case 1:{const_source=false;const_alpha1_value=(T)1.e-3;uvf=false;}break;
    case 2:{const_source=true;const_alpha1_value=(T)1.;uvf=false;}break;
    case 3:{const_source=false;const_alpha1_value=(T)0.;uvf=false;}break;
    case 4:{const_source=true;const_alpha1_value=(T)1.;uvf=false;}break;
    case 5:{const_source=false;const_alpha1_value=(T)0.;uvf=true;}break;
    case 6:{const_source=true;const_alpha1_value=(T)1.;uvf=true;}break;
    case 7:{const_source=true;const_alpha1_value=(T)1.;uvf=false;}break;
    case 8:{const_source=true;const_alpha1_value=(T)1.;uvf=false;}break;}
}
//######################################################################
// Write_Output_Files
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Write_Output_Files(const int frame) const
{
    File_Utilities::Create_Directory(output_directory+"/"+std::to_string(frame));
    File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/frame_title",frame_title);

    // write hierarchy
    File_Utilities::Write_To_Text_File(output_directory+"/"+std::to_string(frame)+"/levels",levels);
    hierarchy->Write_Hierarchy(output_directory,frame);
    hierarchy->template Write_Channel<T>(output_directory+"/"+std::to_string(frame)+"/spgrid_density",alpha1_channel);
    hierarchy->template Write_Channel<T>(output_directory+"/"+std::to_string(frame)+"/spgrid_u",face_velocity1_channels(0));
    hierarchy->template Write_Channel<T>(output_directory+"/"+std::to_string(frame)+"/spgrid_v",face_velocity1_channels(1));
    if(d==3) hierarchy->template Write_Channel<T>(output_directory+"/"+std::to_string(frame)+"/spgrid_w",face_velocity1_channels(2));
}
//######################################################################
// Read_Output_Files
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Read_Output_Files(const int& frame)
{
}
//######################################################################
template class Nova::Spray_Example<float,2>;
template class Nova::Spray_Example<float,3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::Spray_Example<double,2>;
template class Nova::Spray_Example<double,3>;
#endif
