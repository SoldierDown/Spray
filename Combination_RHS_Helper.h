//!#####################################################################
//! \file Combination_RHS_Helper.h
//!#####################################################################
// Class Combination_RHS_Helper
//######################################################################
#ifndef __Combination_RHS_Helper__
#define __Combination_RHS_Helper__


#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
// #include <nova/SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Combination_RHS_Helper
{
    using Hierarchy             = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Combination_RHS_Helper(Hierarchy& hierarchy,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                        T Struct_type::* rhs_channel,T Struct_type::* alpha1_channel,T Struct_type::* alpha2_channel,
                        Channel_Vector& face_velocity1_channels,Channel_Vector& face_velocity2_channels,
                        const T& rho1,const T& rho2,const T& dt,const T& one_over_dx)
    {Run(hierarchy,allocator,blocks,rhs_channel,alpha1_channel,alpha2_channel,face_velocity1_channels,face_velocity2_channels,rho1,rho2,dt,one_over_dx);}

    void Run(Hierarchy& hierarchy,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            T Struct_type::* rhs_channel,T Struct_type::* alpha1_channel,T Struct_type::* alpha2_channel,
            Channel_Vector& face_velocity1_channels,Channel_Vector& face_velocity2_channels,
            const T& rho1,const T& rho2,const T& dt,const T& one_over_dx) const
    {
        const int level=0; const T coeff=(T)1;
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto alpha1=allocator.template Get_Const_Array<Struct_type,T>(alpha1_channel); 
        auto alpha2=allocator.template Get_Const_Array<Struct_type,T>(alpha2_channel); 
        auto rhs=allocator.template Get_Array<Struct_type,T>(rhs_channel);
        uint64_t face_neighbor_offsets[Topology_Helper::number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);
        Vector<uint64_t,d> positive_face_offsets; for(int axis=0;axis<d;++axis) positive_face_offsets(axis)=Topology_Helper::Axis_Vector_Offset(axis);
        Vector<uint64_t,d> negative_face_offsets; for(int axis=0;axis<d;++axis) negative_face_offsets(axis)=Topology_Helper::Negative_Axis_Vector_Offset(axis);
        auto combination_rhs_helper=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) if(flags(offset)&Cell_Type_Interior){T result=(T)0.;
                for(int axis=0;axis<d;++axis){unsigned face_valid_mask=Topology_Helper::Face_Valid_Mask(axis);
                    auto u1=allocator.template Get_Const_Array<Struct_type,T>(face_velocity1_channels(axis));
                    auto u2=allocator.template Get_Const_Array<Struct_type,T>(face_velocity2_channels(axis));
                    // left/down/back face    
                    if(flags(offset)&face_valid_mask) { 
                        const T positive_alpha1=alpha1(offset); const T positive_alpha2=alpha2(offset);
                        if(std::isnan(positive_alpha1)) {Log::cout<<"-NAN positive_alpha1 in Combination_RHS_Helper"<<std::endl;exit(0);}
                        if(std::isnan(positive_alpha2)) {Log::cout<<"-NAN positive_alpha2 in Combination_RHS_Helper"<<std::endl;exit(0);}
                        // interpolate alpha to face
                        uint64_t negative_offset=Flag_array_mask::Packed_Add(offset,negative_face_offsets(axis));
                        T negative_alpha1=(T)0.; T negative_alpha2=(T)0.;
                        if(hierarchy.template Set<unsigned>(level,&Struct_type::flags).Is_Set(negative_offset,Cell_Type_Interior|Cell_Type_Dirichlet)){
                            negative_alpha1=alpha1(negative_offset); negative_alpha2=alpha2(negative_offset);}
                        T interpolated_alpha1=(T).5*(positive_alpha1+negative_alpha1);
                        // T interpolated_alpha2=(T).5*(positive_alpha2+negative_alpha2);
                        T interpolated_alpha2=(T)1.-interpolated_alpha1;
                        // Nova_Utilities::Clamp(interpolated_alpha1,(T)1.e-3,(T)1.); Nova_Utilities::Clamp(interpolated_alpha2,(T)1.e-3,(T)1.);
                        if(std::isnan(u1(offset))) {Log::cout<<"-NAN u1 in Combination_RHS_Helper"<<std::endl;exit(0);}
                        if(std::isnan(u2(offset))) {Log::cout<<"-NAN u2 in Combination_RHS_Helper"<<std::endl;exit(0);}
                        // if(u1(offset)!=u2(offset)) {Log::cout<<"u1 != u2 in Combination_RHS_Helper"<<std::endl;exit(0);}
                        result+=u1(offset)*interpolated_alpha1+u2(offset)*interpolated_alpha2;}
                    
                    // right/up/front face
                    uint64_t positive_offset=Flag_array_mask::Packed_Add(offset,positive_face_offsets(axis));
                    if(flags(positive_offset)&face_valid_mask) {
                        const T positive_alpha1=alpha1(positive_offset); const T positive_alpha2=alpha2(positive_offset);
                        if(std::isnan(positive_alpha1)) {Log::cout<<"+NAN positive_alpha1 in Combination_RHS_Helper"<<std::endl;exit(0);}
                        if(std::isnan(positive_alpha2)) {Log::cout<<"+NAN positive_alpha2 in Combination_RHS_Helper"<<std::endl;exit(0);}
                        // interpolate alpha to face
                        const T negative_alpha1=alpha1(offset); const T negative_alpha2=alpha2(offset);
                        T interpolated_alpha1=(T).5*(positive_alpha1+negative_alpha1);
                        // T interpolated_alpha2=(T).5*(positive_alpha2+negative_alpha2);
                        T interpolated_alpha2=(T)1.-interpolated_alpha1;
                        // Nova_Utilities::Clamp(interpolated_alpha1,(T)1.e-3,(T)1.); Nova_Utilities::Clamp(interpolated_alpha2,(T)1.e-3,(T)1.);
                        if(std::isnan(u1(positive_offset))) {Log::cout<<"+NAN u1 in Combination_RHS_Helper"<<std::endl;exit(0);}
                        if(std::isnan(u2(positive_offset))) {Log::cout<<"+NAN u2 in Combination_RHS_Helper"<<std::endl;exit(0);}
                        // if(u1(positive_offset)!=u2(positive_offset)) {Log::cout<<"u1 != u2 in Combination_RHS_Helper"<<std::endl;exit(0);}
                        result-=u1(positive_offset)*interpolated_alpha1+u2(positive_offset)*interpolated_alpha2;}}
            rhs(offset)=result*one_over_dx*coeff;}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,combination_rhs_helper);
    }

};
}
#endif
