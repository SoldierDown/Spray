//!#####################################################################
//! \file Spray_Example.cpp
//!#####################################################################
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Initializer.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Iterator.h>
#include <nova/Dynamics/Hierarchy/Advection/Grid_Hierarchy_Advection.h>
#include <nova/Tools/Grids/Grid_Iterator_Face.h>
#include <nova/Tools/Krylov_Solvers/Conjugate_Gradient.h>
#include <nova/Tools/Utilities/File_Utilities.h>
#include "Apply_Pressure.h"
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
#include "Relative_Velocity_Helper.h"
#include "Compute_Drag_Force_Helper.h"
#include "Apply_Drag_Force_Helper.h"

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
    rho1=(T)5.;
    rho2=(T)1.;
    face_velocity1_channels(0)                  = &Struct_type::ch0;
    face_velocity1_channels(1)                  = &Struct_type::ch1;
    if(d==3) face_velocity1_channels(2)         = &Struct_type::ch2;
    face_velocity2_channels(0)                  = &Struct_type::ch3;
    face_velocity2_channels(1)                  = &Struct_type::ch4;
    if(d==3) face_velocity2_channels(2)         = &Struct_type::ch5;
    alpha1_channel                              = &Struct_type::ch6;
    alpha1_backup_channel                       = &Struct_type::ch7;
    alpha2_channel                              = &Struct_type::ch8;
    alpha2_backup_channel                       = &Struct_type::ch9;
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
    //Set_Neumann_Faces_Inside_Sources();
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
    Channel_Vector cell_velocity_channels;
    cell_velocity_channels(0)               = &Struct_type::ch10;
    cell_velocity_channels(1)               = &Struct_type::ch11;
    if(d==3) cell_velocity_channels(2)      = &Struct_type::ch12;
    T Struct_type::* temp_channel           = &Struct_type::ch13;
    Vector<uint64_t,d> other_face_offsets;
    for(int axis=0;axis<d;++axis) other_face_offsets(axis)=Topology_Helper::Axis_Vector_Offset(axis);
    Uniform_Grid_Averaging_Helper<Struct_type,T,d>::Uniform_Grid_Average_Face_Velocities_To_Cells(*hierarchy,hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity1_channels,cell_velocity_channels,other_face_offsets);
    Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Density(*hierarchy,cell_velocity_channels,alpha1_channel,temp_channel,dt);
    // ensure alpha_air + alpha_water = 1
    Alpha_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),alpha1_channel,alpha2_channel);
}
//######################################################################
// Backup_Alpha
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Backup_Alpha()
{
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),alpha1_backup_channel);
    SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),alpha2_backup_channel);
    SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),alpha1_channel,
                                        alpha1_backup_channel,(unsigned)Cell_Type_Interior);
    SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),alpha2_channel,
                                        alpha2_backup_channel,(unsigned)Cell_Type_Interior);
}
//######################################################################
// Ficks_Diffusion
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Apply_Drag_Force(const T& dt)
{
    enum {number_of_nodes_per_cell  = Topology_Helper::number_of_nodes_per_cell};
    enum {number_of_nodes_per_face  = Topology_Helper::number_of_nodes_per_face};
    Channel_Vector interpolated_rel_face_velocity_channels;
    interpolated_rel_face_velocity_channels(0)                  = &Struct_type::ch10;
    interpolated_rel_face_velocity_channels(1)                  = &Struct_type::ch11;
    if(d==3) interpolated_rel_face_velocity_channels(2)         = &Struct_type::ch12;
    Channel_Vector rel_face_velocity_channels;
    rel_face_velocity_channels(0)                               = &Struct_type::ch13;
    rel_face_velocity_channels(1)                               = &Struct_type::ch14;
    if(d==3) rel_face_velocity_channels(2)                      = &Struct_type::ch15;
    T Struct_type::* drag_force_channel                         = &Struct_type::ch16;
    Vector<uint64_t,d> negative_face_offsets; for(int axis=0;axis<d;++axis) negative_face_offsets(axis)=Topology_Helper::Negative_Axis_Vector_Offset(axis);
    uint64_t nodes_of_cell_offsets[number_of_nodes_per_cell]; Topology_Helper::Nodes_Of_Cell_Offsets(nodes_of_cell_offsets);
    for(int axis=0;axis<d;++axis) { 
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),rel_face_velocity_channels(axis));
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),drag_force_channel);
        // true: 
        Relative_Velocity_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity1_channels,face_velocity2_channels,rel_face_velocity_channels,true);
        // interpolate velocity at face
        Uniform_Grid_Averaging_Helper<Struct_type,T,d>::Uniform_Grid_Average_Face_Velocities_To_Faces(*hierarchy,hierarchy->Allocator(level),hierarchy->Blocks(level),
                                                    rel_face_velocity_channels,interpolated_rel_face_velocity_channels,negative_face_offsets,nodes_of_cell_offsets,axis);
        Compute_Drag_Force_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),drag_force_channel,alpha1_channel,rel_face_velocity_channels,axis);
        Apply_Drag_Force_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity1_channels(axis),drag_force_channel,alpha1_channel,dt,rho1,axis,true);
    }  
    
}
//######################################################################
// Ficks_Diffusion
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Ficks_Diffusion(const T& dt)
{
    // using Multigrid_struct_type 	= Multigrid_Data<T>;
    // Log::cout<<"Fick's Diffusion"<<std::endl;
	// const Grid<T,d>& grid=hierarchy->Lattice(0);
    // const T one_over_dx2=Nova_Utilities::Sqr(grid.one_over_dX(0));
    // const T a=diff_coeff*dt*one_over_dx2; const T two_d_a_plus_one=(T)2*d*a+(T)1.;
	// if(!explicit_diffusion){
	//     Ficks_CG_System<Struct_type,Multigrid_struct_type,T,d> cg_system(*hierarchy,mg_levels,diff_coeff*dt,3,1,200);
    //     T Struct_type::* q_channel              = &Struct_type::ch8;
    //     T Struct_type::* r_channel              = &Struct_type::ch9;
    //     T Struct_type::* s_channel              = &Struct_type::ch10;
    //     T Struct_type::* k_channel              = &Struct_type::ch11;
    //     T Struct_type::* z_channel              = &Struct_type::ch11;
    //     T Struct_type::* b_channel              = &Struct_type::ch12;
    
    // // clear all channels
    // for(int level=0;level<levels;++level){
    //     SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),q_channel);
    //     SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),r_channel);
    //     SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),s_channel);
    //     SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),z_channel);
    //     SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),b_channel);}

	//     for(int level=0;level<levels;++level) Ficks_RHS_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel,b_channel,a);
	//     for(int level=0;level<levels;++level) SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_backup_channel,
    //                                                                             density_channel,(unsigned)Cell_Type_Interior);
    //     CG_Vector<Struct_type,T,d> x_V(*hierarchy,density_channel),b_V(*hierarchy,b_channel),q_V(*hierarchy,q_channel),
    //                                 s_V(*hierarchy,s_channel),r_V(*hierarchy,r_channel),k_V(*hierarchy,k_channel),z_V(*hierarchy,z_channel);   
	//     Conjugate_Gradient<T> cg;
    //     cg.iterations_used=new int;
    //     cg_system.Multiply(x_V,r_V);
    //     r_V-=b_V;
    //     const T b_norm=cg_system.Convergence_Norm(r_V);
    //     Log::cout<<"Norm: "<<b_norm<<std::endl;
    //     cg.print_residuals=false;
    //     cg.print_diagnostics=true;
    //     cg.restart_iterations=cg_restart_iterations;
    //     const T tolerance=std::max((T)1e-4*b_norm,(T)1e-4);
    //     cg.Solve(cg_system,x_V,b_V,q_V,s_V,r_V,k_V,z_V,tolerance,0,cg_iterations);
    //     iteration_counter+=*(cg.iterations_used);
    //     for(int level=0;level<levels;++level) Density_Clamp_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel);}
    // else{
    //     T Struct_type::* lap_density_channel    = &Struct_type::ch8;
    //     for(int level=0;level<levels;++level) SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),lap_density_channel);
    //     for(int level=0;level<levels;++level) Lap_Calculator<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_backup_channel,lap_density_channel,one_over_dx2);        
    //     for(int level=0;level<levels;++level) SPGrid::Masked_Saxpy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),dt,lap_density_channel,density_backup_channel,density_channel,Cell_Type_Interior);
    //     for(int level=0;level<levels;++level) Density_Clamp_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),density_channel);}
    //     Log::cout<<"Fick's Diffusion finished"<<std::endl;
}
//######################################################################
// Modify_Density_With_Sources
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Modify_Density_With_Sources(){}
//######################################################################
// Add_Source
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Add_Source(const T& dt)
{
    // for(int level=0;level<levels;++level)
    //     Source_Adder<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),density_channel,density_sources,source_rate,dt,level);
}
//######################################################################
// Advect_Face_Velocities
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Advect_Face_Velocities(const T& dt)
{
    Channel_Vector interpolated_face_velocity_channels;
    interpolated_face_velocity_channels(0)              = &Struct_type::ch8;
    interpolated_face_velocity_channels(1)              = &Struct_type::ch9;
    if(d==3) interpolated_face_velocity_channels(2)     = &Struct_type::ch10;
    Channel_Vector face_velocity_backup_channels;
    face_velocity_backup_channels(0)                    = &Struct_type::ch11;
    face_velocity_backup_channels(1)                    = &Struct_type::ch12;
    if(d==3) face_velocity_backup_channels(2)           = &Struct_type::ch13;
    T Struct_type::* temp_channel                       = &Struct_type::ch14;
    // phase 1
    // backup face velocity
    for(int level=0;level<levels;++level) for(int axis=0;axis<d;++axis) {
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity_backup_channels(axis));
        SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity1_channels(axis),
                                                face_velocity_backup_channels(axis),Topology_Helper::Face_Active_Mask(axis));}  
        Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Face_Velocities(*hierarchy,face_velocity1_channels,face_velocity_backup_channels,interpolated_face_velocity_channels,temp_channel,dt);
    // phase 2
    // backup face velocity
    for(int level=0;level<levels;++level) for(int axis=0;axis<d;++axis) {
        SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity_backup_channels(axis));
        SPGrid::Masked_Copy<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),face_velocity2_channels(axis),
                                                face_velocity_backup_channels(axis),Topology_Helper::Face_Active_Mask(axis));}  
        Uniform_Grid_Advection_Helper<Struct_type,T,d>::Uniform_Grid_Advect_Face_Velocities(*hierarchy,face_velocity2_channels,face_velocity_backup_channels,interpolated_face_velocity_channels,temp_channel,dt);
}
//######################################################################
// Set_Neumann_Faces_Inside_Sources
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Set_Neumann_Faces_Inside_Sources()
{
    for(int level=0;level<levels;++level){
        auto flags=hierarchy->Allocator(level).template Get_Array<Struct_type,unsigned>(&Struct_type::flags);
        for(Grid_Iterator_Face<T,d> iterator(hierarchy->Lattice(level));iterator.Valid();iterator.Next()){
            const int axis=iterator.Axis();const T_INDEX& face_index=iterator.Face_Index();
            uint64_t face_offset=Flag_array_mask::Linear_Offset(face_index._data);
            const uint64_t face_active_mask=Topology_Helper::Face_Active_Mask(axis);
            const uint64_t face_valid_mask=Topology_Helper::Face_Valid_Mask(axis);
            const TV X=hierarchy->Lattice(level).Face(axis,face_index);
            for(size_t i=0;i<velocity_sources.size();++i) if(velocity_sources(i)->Inside(X)){
                if(hierarchy->template Set<unsigned>(level,&Struct_type::flags).Is_Set(face_offset,face_valid_mask))
                    flags(face_offset)&=~face_active_mask;
                break;}}}
}
//######################################################################
// Initialize_Velocity_Field
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Initialize_Velocity_Field()
{
    Array<TV> source_velocity;
    source_velocity.Append(TV::Axis_Vector(1)*bv);
    for(int level=0;level<levels;++level)
        if(uvf)Uniform_Velocity_Field_Initializer<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity1_channels,bv,level);
        else Source_Velocity_Setup<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),velocity_sources,source_velocity,face_velocity1_channels,level);
}
//######################################################################
// Project
//######################################################################
template<class T,int d> void Spray_Example<T,d>::
Project()
{
    // using Multigrid_struct_type             = Multigrid_Data<T>;
    // using Hierarchy_Projection              = Grid_Hierarchy_Projection<Struct_type,T,d>;

    // // set up divergence channel
    // T Struct_type::* divergence_channel     = &Struct_type::ch8;
    // T Struct_type::* pressure_channel       = &Struct_type::ch9;

    // // clear
    // for(int level=0;level<levels;++level){
    //     SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),pressure_channel);
    //     SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),divergence_channel);}

    // Array<TV> source_velocity;
    // source_velocity.Append(TV::Axis_Vector(1)*bv);
    // // set boundary conditions
    // for(int level=0;level<levels;++level){Boundary_Condition_Helper<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity_channels,pressure_channel,level);
    //    if(!uvf) Source_Velocity_Setup<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),velocity_sources,source_velocity,face_velocity_channels,level);}
        
    // // compute divergence
    // Hierarchy_Projection::Compute_Divergence(*hierarchy,face_velocity_channels,divergence_channel);

    // Poisson_CG_System<Struct_type,Multigrid_struct_type,T,d> cg_system(*hierarchy,mg_levels,3,1,200);

    // T Struct_type::* q_channel              = &Struct_type::ch10;
    // T Struct_type::* r_channel              = &Struct_type::ch11;
    // T Struct_type::* s_channel              = &Struct_type::ch12;
    // T Struct_type::* k_channel              = &Struct_type::ch12;
    // T Struct_type::* z_channel              = &Struct_type::ch13;

    // // clear all channels
    // for(int level=0;level<levels;++level){
    //     SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),q_channel);
    //     SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),r_channel);
    //     SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),s_channel);
    //     SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),z_channel);}

    // CG_Vector<Struct_type,T,d> x_V(*hierarchy,pressure_channel);
    // CG_Vector<Struct_type,T,d> b_V(*hierarchy,divergence_channel);
    // CG_Vector<Struct_type,T,d> q_V(*hierarchy,q_channel);
    // CG_Vector<Struct_type,T,d> r_V(*hierarchy,r_channel);
    // CG_Vector<Struct_type,T,d> s_V(*hierarchy,s_channel);
    // CG_Vector<Struct_type,T,d> k_V(*hierarchy,k_channel);
    // CG_Vector<Struct_type,T,d> z_V(*hierarchy,z_channel);

    // Conjugate_Gradient<T> cg;
    // cg.print_diagnostics=false; cg.print_residuals=false;
    // cg_system.Multiply(x_V,r_V);
    // r_V-=b_V;
    // const T b_norm=cg_system.Convergence_Norm(r_V);
    // Log::cout<<"Norm: "<<b_norm<<std::endl;
    // cg.restart_iterations=cg_restart_iterations;
    // const T tolerance=std::max((T)1e-4*b_norm,(T)1e-4);
    // cg.Solve(cg_system,x_V,b_V,q_V,s_V,r_V,k_V,z_V,tolerance,0,cg_iterations);

    // for(int level=0;level<levels;++level)
    //     Apply_Pressure<Struct_type,T,d>(*hierarchy,hierarchy->Blocks(level),face_velocity_channels,pressure_channel,level);
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
    parse_args->Add_Double_Argument("-diff_coeff",(T)1e-3,"diffusion coefficient.");
    parse_args->Add_Double_Argument("-fc",(T)0.,"fc.");
    parse_args->Add_Double_Argument("-bv",(T)1.,"Background velocity(along y axis).");
    parse_args->Add_Double_Argument("-sr",(T)1.,"Source rate");
    parse_args->Add_Double_Argument("-tau",(T)1.,"tau.");
    parse_args->Add_Option_Argument("-ficks","Fick's diffusion.");
    parse_args->Add_Option_Argument("-nd","Turn off diffusion");
    parse_args->Add_Option_Argument("-ed","Explicit diffusion");
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
    diff_coeff=parse_args->Get_Double_Value("-diff_coeff");
    bv=parse_args->Get_Double_Value("-bv");
    source_rate=parse_args->Get_Double_Value("-sr");
    Fc=parse_args->Get_Double_Value("-fc");
    tau=parse_args->Get_Double_Value("-tau");
    FICKS=parse_args->Get_Option_Value("-ficks");
    uvf=parse_args->Get_Option_Value("-uvf");
    nd=parse_args->Get_Option_Value("-nd");
    explicit_diffusion=parse_args->Get_Option_Value("-ed");
    cg_iterations=parse_args->Get_Integer_Value("-cg_iterations");
    cg_restart_iterations=parse_args->Get_Integer_Value("-cg_restart_iterations");
    cg_tolerance=(T)parse_args->Get_Double_Value("-cg_tolerance");
    if(nd){diff_coeff=(T)0.;tau=(T)1.;Fc=(T)0.;}
    switch (test_number){
    case 1:{const_source=false;const_alpha1_value=(T)0.;uvf=false;}break;
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
