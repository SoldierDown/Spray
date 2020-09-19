//!#####################################################################
//! \file Standard_Tests.h
//!#####################################################################
// Class Standard_Tests
//######################################################################
#ifndef __Standard_Tests__
#define __Standard_Tests__

#include <nova/Geometry/Implicit_Objects/Box_Implicit_Object.h>
#include <nova/SPGrid/Tools/SPGrid_Clear.h>
#include <nova/Tools/Utilities/Range_Iterator.h>
#include "../../Poisson_Data.h"
#include "../../Spray_Example.h"
#include "../../Rasterizers/Adaptive_Sphere_Rasterizer.h"
#include "../../Rasterizers/Randomized_Rasterizer.h"

#include "../../Clamp_Alpha_Helper.h"

namespace Nova{
template<class T,int d>
class Standard_Tests: public Spray_Example<T,d>
{
    
    using TV                        = Vector<T,d>;
    using T_INDEX                   = Vector<int,d>;
    using Struct_type               = Poisson_Data<T>;
    using Hierarchy                 = Grid_Hierarchy<Struct_type,T,d>;
    using Flags_type                = typename Struct_type::Flags_type;
    using Base                      = Spray_Example<T,d>;
    using Allocator_type            = SPGrid::SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask           = typename Allocator_type::template Array_mask<unsigned>;

  public:
    using Base::output_directory;       using Base::test_number;            using Base::counts; 
    using Base::levels;                 using Base::domain_walls;           using Base::hierarchy;
    using Base::rasterizer;             using Base::cfl;                    using Base::domain;
    using Base::velocity1_sources;      using Base::velocity2_sources;      using Base::use_wall;
    using Base::density_sources;        using Base::inner_wall;
    using Base::rho1;                   using Base::rho2;                   
    using Base::alpha1_channel;         using Base::alpha2_channel;       
    using Base::source_rate;            using Base::uvf;    
    using Base::const_alpha1_value;     using Base::const_source;           
    using Base::explicit_diffusion;
    /****************************
     * example explanation:
     *
     * 1. Simple inflow/outflow.
     ****************************/

    Standard_Tests()
        :Base()
    {}

//######################################################################
    void Parse_Options() override
    {
        Base::Parse_Options();
        output_directory="Spray_2d_rho1_"+std::to_string(rho1)+"_rho2_"+std::to_string(rho2)+"_case_"+std::to_string(test_number)+"_sr_"+std::to_string(source_rate)+"_Resolution_"+std::to_string(counts(0))+"x"+std::to_string(counts(1));
        for(int axis=0;axis<d;++axis) for(int side=0;side<2;++side) domain_walls(axis)(side)=false;
        TV min_corner,max_corner=domain;
        Log::cout<<"domain: "<<domain<<", resolution: "<<counts<<std::endl;
        hierarchy=new Hierarchy(counts,Range<T,d>(min_corner,max_corner),levels);
    }
//######################################################################
    void Initialize_Rasterizer(const int test_number) override
    {
        rasterizer=new Randomized_Rasterizer<Struct_type,T,d>(*hierarchy);
    }
//######################################################################
    void Initialize_Fluid_State(const int test_number) override
    {
        // clear density channel
        for(int level=0;level<levels;++level){
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),alpha1_channel);
            SPGrid::Clear<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),alpha2_channel);
        }

        for(int level=0;level<levels;++level){auto blocks=hierarchy->Blocks(level);
            auto block_size=hierarchy->Allocator(level).Block_Size();
            auto alpha1_data=hierarchy->Allocator(level).template Get_Array<Struct_type,T>(alpha1_channel);
            auto alpha2_data=hierarchy->Allocator(level).template Get_Array<Struct_type,T>(alpha2_channel);
            auto flags=hierarchy->Allocator(level).template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);

            for(unsigned block=0;block<blocks.second;++block){uint64_t offset=blocks.first[block];
                Range_Iterator<d> range_iterator(T_INDEX(),*reinterpret_cast<T_INDEX*>(&block_size)-1);
                T_INDEX base_index(Flag_array_mask::LinearToCoord(offset));

                for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)){
                    const T_INDEX index=base_index+range_iterator.Index();
                    if(flags(offset)&Cell_Type_Interior) {alpha1_data(offset)=const_alpha1_value; alpha2_data(offset)=(T)1.-const_alpha1_value;}
                    range_iterator.Next();}}}
        const int level=0;
        // Clamp_Alpha_Helper<Struct_type,T,d>(hierarchy->Allocator(level),hierarchy->Blocks(level),alpha1_channel,alpha2_channel);
        
    }
//######################################################################
    void Initialize_Sources(const int test_number) override
    {
        const T cell_width=domain(0)/counts(0);
        switch (test_number)
        {
        // test case 1: density&velocity source near the bottom 
        case 1:
        case 2:{
            TV density_min_corner=TV({(T)1.8,(T)0.}),density_max_corner=TV({(T)2.2,(T)2.*cell_width});
            Implicit_Object<T,d>* density_obj=new Box_Implicit_Object<T,d>(density_min_corner,density_max_corner);
            density_sources.Append(density_obj);

            TV velocity_min_corner=TV({(T)1.8,(T)0.}),velocity_max_corner=TV({(T)2.2,(T)2.*cell_width});
            Implicit_Object<T,d>* velocity_obj=new Box_Implicit_Object<T,d>(velocity_min_corner,velocity_max_corner);
            velocity1_sources.Append(velocity_obj);
            
            if(use_wall){}
            }break;
        case 3:
        case 4:{
            TV density_min_corner=TV({(T)2.-cell_width,(T)2.-cell_width}),density_max_corner=TV({(T)2.+cell_width,(T)2.+cell_width});
            Implicit_Object<T,d>* density_obj=new Box_Implicit_Object<T,d>(density_min_corner,density_max_corner);
            density_sources.Append(density_obj);
            TV velocity_min_corner=TV({(T)1.8,(T)0.}),velocity_max_corner=TV({(T)2.2,(T)2.*cell_width});
            Implicit_Object<T,d>* velocity_obj=new Box_Implicit_Object<T,d>(velocity_min_corner,velocity_max_corner);
            velocity1_sources.Append(velocity_obj);
        }break;
        case 5:
        case 6:
        case 7:
        case 8:{
            TV density_min_corner=TV({(T)2.-cell_width,(T)2.-cell_width}),density_max_corner=TV({(T)2.+cell_width,(T)2.+cell_width});
            Implicit_Object<T,d>* density_obj=new Box_Implicit_Object<T,d>(density_min_corner,density_max_corner);
            density_sources.Append(density_obj);
            TV velocity_min_corner=TV({(T)2.-cell_width,(T)2.-cell_width}),velocity_max_corner=TV({(T)2.+cell_width,(T)2.+cell_width});
            Implicit_Object<T,d>* velocity_obj=new Box_Implicit_Object<T,d>(velocity_min_corner,velocity_max_corner);
            velocity1_sources.Append(velocity_obj);}break;
        case 9:{
            TV density_min_corner=TV({(T)6.4,(T)4.-(T)2.*cell_width}),density_max_corner=TV({(T)6.8,(T)4.});
            Implicit_Object<T,d>* density_obj=new Box_Implicit_Object<T,d>(density_min_corner,density_max_corner);
            density_sources.Append(density_obj);

            TV velocity_min_corner=TV({(T)8.-(T)2.*cell_width,(T)0.}),velocity_max_corner=TV({(T)8.,(T)4.});
            Implicit_Object<T,d>* velocity_obj=new Box_Implicit_Object<T,d>(velocity_min_corner,velocity_max_corner);
            velocity1_sources.Append(velocity_obj);
        }break;}

    }
//######################################################################
};
}
#endif
