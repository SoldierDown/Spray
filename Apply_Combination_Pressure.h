//!#####################################################################
//! \file Apply_Combination_Pressure.h
//!#####################################################################
// Class Apply_Combination_Pressure
//######################################################################
#ifndef __Apply_Combination_Pressure__
#define __Apply_Combination_Pressure__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Apply_Combination_Pressure
{
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Apply_Combination_Pressure(Hierarchy& hierarchy,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                            Channel_Vector& face_velocity_channels,T Struct_type::* pressure_channel,const T& rho,const T& one_over_dx)
    {Run(hierarchy,allocator,blocks,face_velocity_channels,pressure_channel,rho,one_over_dx);}

    void Run(Hierarchy& hierarchy,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            Channel_Vector& face_velocity_channels,T Struct_type::* pressure_channel,const T& rho,const T& one_over_dx) const
    {
        const int level=0;
        auto pressure=allocator.template Get_Const_Array<Struct_type,T>(pressure_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        Vector<uint64_t,d> negative_face_offsets;
        for(int axis=0;axis<d;++axis) negative_face_offsets(axis)=Topology_Helper::Negative_Axis_Vector_Offset(axis);

        auto apply_combination_pressure=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) for(int axis=0;axis<d;++axis){
                auto face_velocity=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(face_velocity_channels(axis));
                if(flags(offset)&Topology_Helper::Face_Active_Mask(axis)){uint64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,negative_face_offsets(axis));
                    T one_over_dX=hierarchy.Lattice(level).one_over_dX[axis];
                    face_velocity(offset)-=(pressure(offset)-pressure(neighbor_offset))*one_over_dX*rho;}}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,apply_combination_pressure);
    }
};
}
#endif
