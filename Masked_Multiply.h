//!#####################################################################
//! \file Masked_Multiply.h
//!#####################################################################
// Class Masked_Multiply
//######################################################################
#ifndef __Masked_Multiply__
#define __Masked_Multiply__

#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>

namespace SPGrid{
template<class Struct_type,class T,int d>
class Masked_Multiply
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Allocator_type        = SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;

  public:
    Masked_Multiply(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* source1_channel,
               T Struct_type::* source2_channel,T Struct_type::* destination_channel,const unsigned mask)
    {Run(allocator,blocks,source1_channel,source2_channel,destination_channel,mask);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* source1_channel,
             T Struct_type::* source2_channel,T Struct_type::* destination_channel,const unsigned mask) const
    {
        auto source1=allocator.template Get_Const_Array<Struct_type,T>(source1_channel);
        auto source2=allocator.template Get_Const_Array<Struct_type,T>(source2_channel);
        auto destination=allocator.template Get_Array<Struct_type,T>(destination_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        auto masked_multiply=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type))
                if(flags(offset)&mask) destination(offset)=source1(offset)*source2(offset);
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,masked_multiply);
    }
};
}
#endif