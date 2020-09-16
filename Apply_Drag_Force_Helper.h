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
        T Struct_type::* velocity_channel, T Struct_type::* drag_force_channel, T Struct_type::* alpha_channel,
        const T& dt, const T& rho, const int& axis,const bool& f12)
    { Run(hierarchy,blocks,velocity_channel,drag_force_channel,alpha_channel,dt,rho,axis,f12);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
        T Struct_type::* velocity_channel, T Struct_type::* drag_force_channel, T Struct_type::* alpha_channel,
        const T& dt, const T& rho, const int& axis,const bool& f12) const
    {
        const int level=0; const T sign=f12?(T)1.:(T)-1.;
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto v_data=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(velocity_channel);
        auto force_data=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(drag_force_channel);
        auto alpha_data=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(alpha_channel);
        

        auto apply_drag_force_helper=[&](uint64_t offset)
        {
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                    const unsigned face_valid_mask=Topology_Helper::Face_Valid_Mask(axis); const unsigned face_active_mask=Topology_Helper::Face_Active_Mask(axis);
                    if(flags(offset)&face_valid_mask)  v_data(offset)+=sign*dt/(alpha_data(offset)*rho)*force_data(offset);                    
                range_iterator.Next();}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,apply_drag_force_helper);
    }

};
}
#endif
