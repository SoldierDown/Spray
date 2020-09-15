//!#####################################################################
//! \file Spray_Driver.h
//!#####################################################################
// Class Spray_Driver
//######################################################################
#ifndef __Spray_Driver__
#define __Spray_Driver__

#include <nova/Tools/Utilities/Driver.h>
#include "Spray_Example.h"

namespace Nova{
template<class T,int d>
class Spray_Driver: public Driver<T,d>
{
    using TV                = Vector<T,d>;
    using Base              = Driver<T,d>;

    using Base::time;
    using Base::Compute_Dt;using Base::Write_Output_Files;

  public:
    int substep_counter;
    T density_advection_rt;
    T velocity_advection_rt;
    T source_modification_rf;
    T projection_rt;
    T total_rt;
    Spray_Example<T,d>& example;

    Spray_Driver(Spray_Example<T,d>& example_input);
    ~Spray_Driver() {}

//######################################################################
    void Initialize() override;
    void Advance_One_Time_Step_Explicit_Part(const T dt,const T time);
    void Advance_One_Time_Step_Implicit_Part(const T dt,const T time);
    void Advance_To_Target_Time(const T target_time) override;
    void Simulate_To_Frame(const int frame) override;
//######################################################################
};
}
#endif
