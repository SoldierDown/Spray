//!#####################################################################
//! \file Combination_Multiply_Helper.h
//!#####################################################################
// Class Combination_Multiply_Helper
//######################################################################
#ifndef __Combination_Multiply_Helper__
#define __Combination_Multiply_Helper__

#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
#include <nova/Dynamics/Hierarchy/Grid_Topology_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Tools/Utilities/Range_Iterator.h>

namespace Nova{
template<class Struct_type,class T,int d>
class Combination_Multiply_Helper
{
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Channel_Vector            = Vector<T Struct_type::*,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;
    using Topology_Helper           = Grid_Topology_Helper<Flag_array_mask>;
    enum {number_of_faces_per_cell  = Topology_Helper::number_of_faces_per_cell};

  public:
    Combination_Multiply_Helper(Hierarchy& hierarchy,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
                            T Struct_type::* p_channel,T Struct_type::* alpha_channel,T Struct_type::* Ap_channel,
                            T rho1,T rho2)
    {Run(hierarchy,allocator,blocks,p_channel,alpha_channel,Ap_channel,rho1,rho2);}

    void Run(Hierarchy& hierarchy,Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,
            T Struct_type::* p_channel,T Struct_type::* alpha_channel,T Struct_type::* Ap_channel,
            T rho1,T rho2) const
    {
        const int level=0;
        auto block_size=allocator.Block_Size();
        auto pressure=allocator.template Get_Const_Array<Struct_type,T>(p_channel);
        auto alpha=allocator.template Get_Const_Array<Struct_type,T>(&Struct_type::ch12);
        auto Ap_data=allocator.template Get_Array<Struct_type,T>(Ap_channel);
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

        const T c=rho1-rho2;
        uint64_t face_neighbor_offsets[number_of_faces_per_cell];
        Topology_Helper::Face_Neighbor_Offsets(face_neighbor_offsets);

        Vector<uint64_t,d> positive_face_offsets; for(int axis=0;axis<d;++axis) positive_face_offsets(axis)=Topology_Helper::Axis_Vector_Offset(axis);
        Vector<uint64_t,d> negative_face_offsets; for(int axis=0;axis<d;++axis) negative_face_offsets(axis)=Topology_Helper::Negative_Axis_Vector_Offset(axis);

        const T one_over_dx2=Nova_Utilities::Sqr<T>(hierarchy.Lattice(level).one_over_dX[0]);
        const T one_over_2dx=(T).5*hierarchy.Lattice(level).one_over_dX(0);
        auto combination_multiply_helper=[&](uint64_t offset)
        {
            for(unsigned e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                if(flags(offset)&Cell_Type_Interior){
                    const T coeff1=rho1-c*alpha(offset); const T coeff2=c; T laplacian=(T)0.; T gradient=(T)0.;
                    // Log::cout<<"coeff1: "<<coeff1<<", coeff2: "<<coeff2<<std::endl;
                    for(int face=0;face<number_of_faces_per_cell;++face){
                        uint64_t neighbor_offset=Flag_array_mask::Packed_Add(offset,face_neighbor_offsets[face]);
                        if(hierarchy.template Set<unsigned>(level,&Struct_type::flags).Is_Set(neighbor_offset,Cell_Type_Interior)) 
                            laplacian+=one_over_dx2*(pressure(offset)-pressure(neighbor_offset));
                        else if(hierarchy.template Set<unsigned>(level,&Struct_type::flags).Is_Set(neighbor_offset,Cell_Type_Dirichlet)) 
                            laplacian+=one_over_dx2*pressure(offset);} 
                    if(std::isnan(laplacian)) {Log::cout<<"NAN laplacian in Combination_Multiply_Helper"<<std::endl;exit(0);}
                    laplacian*=coeff1;
                    for(int axis=0;axis<d;++axis) { T gradient_pressure=(T)0.; T gradient_alpha=(T)0.;
                        uint64_t positive_offset=Flag_array_mask::Packed_Add(offset,positive_face_offsets(axis));
                        if(hierarchy.template Set<unsigned>(level,&Struct_type::flags).Is_Set(positive_offset,Cell_Type_Interior)) {
                            gradient_pressure+=one_over_2dx*pressure(positive_offset);
                            gradient_alpha+=one_over_2dx*alpha(positive_offset);
                            if(std::isnan(gradient_pressure)) {Log::cout<<"+NAN gradient_pressure in Combination_Multiply_Helper"<<std::endl;exit(0);}
                            if(std::isnan(gradient_alpha)) {Log::cout<<"+NAN gradient_alpha in Combination_Multiply_Helper"<<std::endl;exit(0);}
                
                            }
                        uint64_t negative_offset=Flag_array_mask::Packed_Add(offset,negative_face_offsets(axis));
                        if(hierarchy.template Set<unsigned>(level,&Struct_type::flags).Is_Set(negative_offset,Cell_Type_Interior)) {
                            gradient_pressure+=-one_over_2dx*pressure(negative_offset);
                            gradient_alpha+=-one_over_2dx*alpha(negative_offset);
                            if(std::isnan(gradient_pressure)) {Log::cout<<"-NAN gradient_pressure in Combination_Multiply_Helper"<<std::endl;exit(0);}
                            if(std::isnan(gradient_alpha)) {Log::cout<<"-NAN gradient_alpha in Combination_Multiply_Helper"<<std::endl;exit(0);}
                            }
                        
                        gradient+=gradient_pressure*gradient_alpha;}
                    gradient*=coeff2;
                    // Log::cout<<"laplacian: "<<laplacian<<", gradient: "<<gradient<<std::endl;

                    Ap_data(offset)=laplacian+gradient;}}
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,combination_multiply_helper);
    }
};
}
#endif
