//!#####################################################################
//! \file Exp_Helper.h
//!#####################################################################
// Class Exp_Helper
//######################################################################
#ifndef __Exp_Helper__
#define __Exp_Helper__

#include <nova/Geometry/Implicit_Objects/Implicit_Object.h>
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Exp_Helper
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
    Exp_Helper(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
        T Struct_type::* velocity1_channel,T Struct_type::* velocity2_channel,T Struct_type::* alpha1_channel,
        const T& dt, const T& rho1,const T& rho2,const int& axis)
    {Run(hierarchy,blocks,velocity1_channel,velocity2_channel,alpha1_channel,dt,rho1,rho2,axis);}

    void Run(Hierarchy& hierarchy,const std::pair<const uint64_t*,unsigned>& blocks,
        T Struct_type::* velocity1_channel,T Struct_type::* velocity2_channel,T Struct_type::* alpha1_channel,
        const T& dt, const T& rho1,const T& rho2,const int& axis) const
    {
        const int level=0; 
        auto block_size=hierarchy.Allocator(level).Block_Size();
        auto flags=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto alpha1_data=hierarchy.Allocator(level).template Get_Const_Array<Struct_type,T>(alpha1_channel);
        auto v1_data=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(velocity1_channel);
        auto v2_data=hierarchy.Allocator(level).template Get_Array<Struct_type,T>(velocity2_channel);
        
        const T mu=(T)1.; const T r=1.e-2;

        Vector<uint64_t,d> negative_face_offsets;
        for(int axis=0;axis<d;++axis) negative_face_offsets(axis)=Topology_Helper::Negative_Axis_Vector_Offset(axis);
        
        auto exp_helper=[&](uint64_t offset)
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
                        const T beta=(T)9./(T)2.*interpolated_alpha1*rho2*mu/Nova_Utilities::Sqr(r);
                        const T beta1=beta/(interpolated_alpha1*rho1); const T beta2=beta/((interpolated_alpha2*rho2));
                        const T coeff1=std::exp(-beta1*dt); const T coeff2=std::exp(-beta2*dt);
                        const T u1=v1_data(offset); const T u2=v2_data(offset);

                        v2_data(offset)=u1+coeff2*(u2-u1);
                        v1_data(offset)=u2-coeff1*(u2-u1);

                        if(std::isnan(interpolated_alpha1)) {Log::cout<<"NAN interpolated_alpha1 in Apply_Drag_Force_Helper"<<std::endl;exit(0);}
                        if(std::isnan(interpolated_alpha2)) {Log::cout<<"NAN interpolated_alpha2 in Apply_Drag_Force_Helper"<<std::endl;exit(0);}
                        if(std::isnan(v1_data(offset))) {Log::cout<<"NAN v1 in Apply_Drag_Force_Helper"<<std::endl;exit(0);}
                        if(std::isnan(v2_data(offset))) {Log::cout<<"NAN v2 in Apply_Drag_Force_Helper"<<std::endl;exit(0);}
                }            
                range_iterator.Next();}
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,exp_helper);
    }

};
}
#endif
