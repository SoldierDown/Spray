//!#####################################################################
//! \file Apply_External_Force_Helper.h
//!#####################################################################
// Class Apply_External_Force_Helper
//######################################################################
#ifndef __Apply_External_Force_Helper__
#define __Apply_External_Force_Helper__

#include <nova/Geometry/Implicit_Objects/Implicit_Object.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Apply_External_Force_Helper
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
    Apply_External_Force_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Channel_Vector& face_velocity_channels,const TV& gravity,const T& dt)
    {
        Run(allocator,blocks,face_velocity_channels,gravity,dt);
    }

    void Run(SPGrid::SPGrid_Allocator<Struct_type,2>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Vector<T Struct_type::*,2>& face_velocity_channels,const Vector<T,2>& gravity,const T& dt) const
    {
        auto block_size=allocator.Block_Size();
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto vx=allocator.template Get_Array<Struct_type,T>(face_velocity_channels(0));
        auto vy=allocator.template Get_Array<Struct_type,T>(face_velocity_channels(1));

        auto apply_external_force_helper=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) for(int axis=0;axis<d;++axis){
                auto face_velocity=allocator.template Get_Array<Struct_type,T>(face_velocity_channels(axis));
                if(flags(offset)&Topology_Helper::Face_Active_Mask(axis)){
                    vx(offset)+=dt*gravity(0); vy(offset)+=dt*gravity(1);
                    if(std::isnan(vx(offset))) {Log::cout<<"NAN vx in Apply_External_Force_Helper"<<std::endl;exit(0);}
                    if(std::isnan(vy(offset))) {Log::cout<<"NAN vy in Apply_External_Force_Helper"<<std::endl;exit(0);}
                    
                    }}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,apply_external_force_helper);
    }

    void Run(SPGrid::SPGrid_Allocator<Struct_type,3>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,Vector<T Struct_type::*,3>& face_velocity_channels,const Vector<T,3>& gravity,const T& dt) const
    {
        auto block_size=allocator.Block_Size();
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto vx=allocator.template Get_Array<Struct_type,T>(face_velocity_channels(0));
        auto vy=allocator.template Get_Array<Struct_type,T>(face_velocity_channels(1));
        auto vz=allocator.template Get_Array<Struct_type,T>(face_velocity_channels(2));

        auto apply_external_force_helper=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) for(int axis=0;axis<d;++axis){
                auto face_velocity=allocator.template Get_Array<Struct_type,T>(face_velocity_channels(axis));
                if(flags(offset)&Topology_Helper::Face_Active_Mask(axis)){
                    vx(offset)+=dt*gravity(0); vy(offset)+=dt*gravity(1); vz(offset)+=dt*gravity(2);}}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,apply_external_force_helper);
    }

};
}
#endif
