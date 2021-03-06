//!#####################################################################
//! \file Apply_Drag_Force_Helper.h
//!#####################################################################
// Class Apply_Drag_Force_Helper
//######################################################################
#ifndef __Apply_Drag_Force_Helper__
#define __Apply_Drag_Force_Helper__

#include <nova/Geometry/Implicit_Objects/Implicit_Object.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Apply_Drag_Force_Helper
{
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Apply_Drag_Force_Helper(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
        T Struct_type::* velocity1_channel,T Struct_type::* velocity2_channel, T Struct_type::* drag_force_channel, 
        T Struct_type::* alpha1_channel,T Struct_type::* alpha2_channel,const T& dt, const T& rho1,const T& rho2, const int& axis)
    {Run(hierarchy,blocks,velocity1_channel,velocity2_channel,drag_force_channel,alpha1_channel,alpha2_channel,dt,rho1,rho2,axis);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
        T Struct_type::* velocity1_channel,T Struct_type::* velocity2_channel, T Struct_type::* drag_force_channel, 
        T Struct_type::* alpha1_channel,T Struct_type::* alpha2_channel,const T& dt,const T& rho1,const T& rho2,const int& axis) const
    {
        const int level=0; 
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto force_data=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(drag_force_channel);
        auto alpha1_data=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(alpha1_channel);
        auto alpha2_data=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(alpha2_channel);
        auto v1_data=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(velocity1_channel);
        auto v2_data=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(velocity2_channel);
        
        Vector<uint64_t,d> negative_face_offsets;
        for(int axis=0;axis<d;++axis) negative_face_offsets(axis)=Topology_Helper::Negative_Axis_Vector_Offset(axis);
        
        auto apply_drag_force_helper=[&](uint64_t offset)
        {
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                    const unsigned face_active_mask=Topology_Helper::Face_Active_Mask(axis);
                    if(flags(offset)&face_active_mask) { 
                        uint64_t negative_offset=Flag_array_mask::Packed_Add(offset,negative_face_offsets(axis));
                        T interpolated_alpha1=(T).5*(alpha1_data(offset)+alpha1_data(negative_offset));
                        interpolated_alpha1=Nova_Utilities::Clamp(interpolated_alpha1,(T)1.e-3,(T)1.-(T)1.e-3);
                        T interpolated_alpha2=(T)1.-interpolated_alpha1;
                        v1_data(offset)+=dt/(interpolated_alpha1*rho1)*force_data(offset);
                        v2_data(offset)=dt/(interpolated_alpha2*rho2)*force_data(offset);
                        if(std::isnan(interpolated_alpha1)) {Log::cout<<"NAN interpolated_alpha1 in Apply_Drag_Force_Helper"<<std::endl;exit(0);}
                        if(std::isnan(interpolated_alpha2)) {Log::cout<<"NAN interpolated_alpha2 in Apply_Drag_Force_Helper"<<std::endl;exit(0);}
                        if(std::isnan(force_data(offset))) {Log::cout<<"NAN force in Apply_Drag_Force_Helper"<<std::endl;exit(0);}
                        if(std::isnan(v1_data(offset))) {Log::cout<<"NAN v1 in Apply_Drag_Force_Helper"<<std::endl;exit(0);}
                        if(std::isnan(v2_data(offset))) {Log::cout<<"NAN v2 in Apply_Drag_Force_Helper"<<std::endl;exit(0);}
                }            
                range_iterator.Next();}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,apply_drag_force_helper);
    }

};
}
#endif
