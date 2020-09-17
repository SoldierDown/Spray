//!#####################################################################
//! \file Clamp_Alpha_Helper.h
//!#####################################################################
// Class Clamp_Alpha_Helper
//######################################################################
#ifndef __Clamp_Alpha_Helper__
#define __Clamp_Alpha_Helper__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Tools/Vectors/Vector.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Clamp_Alpha_Helper
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Channel_Vector        = Vector<T Struct_type::*,d>;
    using Allocator_type        = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper       = Grid_Topology_Helper<Flag_array_mask>;

  public:
    Clamp_Alpha_Helper(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                T Struct_type::* alpha1_channel,T Struct_type::* alpha2_channel)
    {Run(allocator,blocks,alpha1_channel,alpha2_channel);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                T Struct_type::* alpha1_channel,T Struct_type::* alpha2_channel) const
    {
        Log::cout<<"clamping ..."<<std::endl;
        auto alpha1_data=allocator.template Get_Array<Struct_type,T>(alpha1_channel);
        auto alpha2_data=allocator.template Get_Array<Struct_type,T>(alpha2_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto clamp_alpha_helper=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&Cell_Type_Interior) { 
                    const T threshold=(T)1.e-3; const T alpha1=alpha1_data(offset); const T alpha2=alpha2_data(offset);
                    // if((alpha1+alpha2)!=(T)1.) Log::cout<<"alpha1 + alpha2 != 1: "<<alpha1<<","<<alpha2<<std::endl;
                    if(alpha1<threshold){alpha1_data(offset)=threshold;alpha2_data(offset)=(T)1.-threshold;}
                    else if(alpha2<threshold){alpha2_data(offset)=threshold; alpha1_data(offset)=(T)1.-threshold;}
                }
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,clamp_alpha_helper);
    }

};
}
#endif