//!#####################################################################
//! \file Check_Velocity_Helper.h
//!#####################################################################
// Class Check_Velocity_Helper
//######################################################################
#ifndef __Check_Velocity_Helper__
#define __Check_Velocity_Helper__

#include <nova/Geometry/Implicit_Objects/Implicit_Object.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Check_Velocity_Helper
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
    Check_Velocity_Helper(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
                        Channel_Vector& face_velocity1_channels,Channel_Vector& face_velocity2_channels)
    {Run(hierarchy,blocks,face_velocity1_channels,face_velocity2_channels);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
            Channel_Vector& face_velocity1_channels,Channel_Vector& face_velocity2_channels) const
    {
        const int level=0;
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto check_velocity_helper=[&](uint64_t offset)
        {
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                for(int axis=0;axis<d;++axis) if(flags(offset)&Topology_Helper::Face_Valid_Mask(axis)){
                    auto v1=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(face_velocity1_channels(axis));
                    auto v2=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(face_velocity2_channels(axis));
                    if(std::abs(v1(offset)-v2(offset))>1e-5*std::abs(v1(offset))){Log::cout<<"v1 != v2 in Check_Velocity_Helper: "<<v1(offset)<<","<<v2(offset)<<std::endl;exit(0);}
                    
                }

                range_iterator.Next();}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,check_velocity_helper);
    }
};
}
#endif
