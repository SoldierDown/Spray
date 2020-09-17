//!#####################################################################
//! \file Relative_Velocity_Helper.h
//!#####################################################################
// Class Relative_Velocity_Helper
//######################################################################
#ifndef __Relative_Velocity_Helper__
#define __Relative_Velocity_Helper__

#include <nova/Geometry/Implicit_Objects/Implicit_Object.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Relative_Velocity_Helper
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
    Relative_Velocity_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
        Channel_Vector& face_velocity1_channels,Channel_Vector& face_velocity2_channels, Channel_Vector& rel_face_velocity_channels,
        const bool& u12)
    { Run(allocator,blocks,face_velocity1_channels,face_velocity2_channels,rel_face_velocity_channels,u12);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
        Channel_Vector& face_velocity1_channels,Channel_Vector& face_velocity2_channels, Channel_Vector& rel_face_velocity_channels,
        const bool& u12) const
    {
        const int level=0;
        auto block_size=allocator.Block_Size();
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        auto u12_helper=[&](uint64_t offset)
        {
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                for(int axis=0;axis<d;++axis){
                    auto v1_data=allocator.template Get_Const_Array<Struct_type,T>(face_velocity1_channels(axis));
                    auto v2_data=allocator.template Get_Const_Array<Struct_type,T>(face_velocity2_channels(axis));
                    auto v12_data=allocator.template Get_Array<Struct_type,T>(rel_face_velocity_channels(axis));
                    const unsigned face_valid_mask=Topology_Helper::Face_Valid_Mask(axis); const unsigned face_active_mask=Topology_Helper::Face_Active_Mask(axis);
                    if(flags(offset)&face_valid_mask)  {
                        v12_data(offset)=v2_data(offset)-v1_data(offset);
                        if(std::isnan(v2_data(offset))) {Log::cout<<"NAN v1 in Relative_Velocity_Helper"<<std::endl;exit(0);}
                        if(std::isnan(v1_data(offset))) {Log::cout<<"NAN v2 in Relative_Velocity_Helper"<<std::endl;exit(0);}
                    }
                    
                    
                    }
                range_iterator.Next();}
        };
        auto u21_helper=[&](uint64_t offset)
        {
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                for(int axis=0;axis<d;++axis){
                    auto v1_data=allocator.template Get_Const_Array<Struct_type,T>(face_velocity1_channels(axis));
                    auto v2_data=allocator.template Get_Const_Array<Struct_type,T>(face_velocity2_channels(axis));
                    auto v21_data=allocator.template Get_Array<Struct_type,T>(rel_face_velocity_channels(axis));
                    const unsigned face_valid_mask=Topology_Helper::Face_Valid_Mask(axis); const unsigned face_active_mask=Topology_Helper::Face_Active_Mask(axis);
                    if(flags(offset)&face_valid_mask)  v21_data(offset)=v1_data(offset)-v2_data(offset);}
                range_iterator.Next();}
        };
        if(u12) SPGrid_Computations::Run_Parallel_Blocks(blocks,u12_helper);
        else SPGrid_Computations::Run_Parallel_Blocks(blocks,u21_helper);
    }
};
}
#endif
