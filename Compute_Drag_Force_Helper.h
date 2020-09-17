//!#####################################################################
//! \file Compute_Drag_Force_Helper.h
//!#####################################################################
// Class Compute_Drag_Force_Helper
//######################################################################
#ifndef __Compute_Drag_Force_Helper__
#define __Compute_Drag_Force_Helper__

#include <nova/Geometry/Implicit_Objects/Implicit_Object.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Compute_Drag_Force_Helper
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
    Compute_Drag_Force_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
        T Struct_type::* drag_force_channel, T Struct_type::* alpha_channel, Channel_Vector& rel_face_velocity_channels,
        const T& rho2,const int& axis)
    { Run(allocator,blocks,drag_force_channel,alpha_channel,rel_face_velocity_channels,rho2,axis);}

    void Run(SPGrid::SPGrid_Allocator<Struct_type,2>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
        T Struct_type::* drag_force_channel, T Struct_type::* alpha_channel, Vector<T Struct_type::*,2>& rel_face_velocity_channels,
        const T& rho2,const int& axis) const
    {
        const int level=0;
        const T D=(T)1.e-2;
        const T mu=(T)1.;
        auto block_size=allocator.Block_Size();
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto vx_data=allocator.template Get_Const_Array<Struct_type,T>(rel_face_velocity_channels(0));
        auto vy_data=allocator.template Get_Const_Array<Struct_type,T>(rel_face_velocity_channels(1));
        auto v_axis_data=allocator.template Get_Const_Array<Struct_type,T>(rel_face_velocity_channels(axis));
        auto alpha_data=allocator.template Get_Const_Array<Struct_type,T>(alpha_channel);        
        auto force_data=allocator.template Get_Array<Struct_type,T>(drag_force_channel);

        auto compute_drag_force_helper=[&](uint64_t offset)
        {
            Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
            T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                    const unsigned face_active_mask=Topology_Helper::Face_Active_Mask(axis);
                    if(flags(offset)&face_active_mask) { 
                        T v_axis_value=v_axis_data(offset); T v_norm=std::sqrt(Nova_Utilities::Sqr(vx_data(offset))+Nova_Utilities::Sqr(vy_data(offset)));
                        force_data(offset)=(T)9./(T)2.*alpha_data(offset)*rho2*mu/Nova_Utilities::Sqr(D)*v_axis_value;
                        if(std::isnan(v_axis_value)) {Log::cout<<"NAN v_axis_value in Compute_Drag_Force_Helper"<<std::endl;exit(0);}
                        if(std::isnan(v_norm)) {Log::cout<<"NAN v_norm in Compute_Drag_Force_Helper"<<std::endl;exit(0);}
                        if(std::isnan(alpha_data(offset))) {Log::cout<<"NAN alpha in Compute_Drag_Force_Helper"<<std::endl;exit(0);}
                        
                        
                        }
                range_iterator.Next();}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,compute_drag_force_helper);
    }
    void Run(SPGrid::SPGrid_Allocator<Struct_type,3>& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
        T Struct_type::* drag_force_channel, T Struct_type::* alpha_channel, Vector<T Struct_type::*,3>& rel_face_velocity_channels,
        const T& rho2,const int& axis) const
    {
        // const int level=0;
        // const T D=(T)1.e-2;
        // const T rho=(T)5.;
        // const T mu=(T)1.;
        // auto block_size=allocator.Block_Size();
        // auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        // auto vx_data=allocator.template Get_Const_Array<Struct_type,T>(rel_face_velocity_channels(0));
        // auto vy_data=allocator.template Get_Const_Array<Struct_type,T>(rel_face_velocity_channels(1));
        // auto vz_data=allocator.template Get_Const_Array<Struct_type,T>(rel_face_velocity_channels(2));
        // auto v_axis_data=allocator.template Get_Array<Struct_type,T>(rel_face_velocity_channels(axis));
        // auto force_data=allocator.template Get_Array<Struct_type,T>(drag_force_channel);
        // auto alpha_data=allocator.template Get_Const_Array<Struct_type,T>(alpha_channel);
        

        // auto compute_drag_force_helper=[&](uint64_t offset)
        // {
        //     Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
        //     T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));
        //     for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
        //             const unsigned face_valid_mask=Topology_Helper::Face_Valid_Mask(axis); const unsigned face_active_mask=Topology_Helper::Face_Active_Mask(axis);
        //             if(flags(offset)&face_valid_mask) { T v_axis_value=v_axis_data(offset); 
        //                 T v_norm=std::sqrt(Nova_Utilities::Sqr(vx_data(offset))+Nova_Utilities::Sqr(vy_data(offset))+Nova_Utilities::Sqr(vz_data(offset)));
        //                 force_data(offset)=(T)9.*std::pow((T).5,gamma)*alpha_data(offset)*rho*std::pow(mu,gamma)*std::pow(D,-(1+gamma))*std::pow(v_norm,1-gamma)*v_axis_value;}
        //         range_iterator.Next();}
        // };
        // SPGrid_Computations::Run_Parallel_Blocks(blocks,compute_drag_force_helper);
    }
};
}
#endif
