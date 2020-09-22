//!#####################################################################
//! \file Poisson_CG_System.h
//!#####################################################################
// Class Poisson_CG_System
//######################################################################
#ifndef __Poisson_CG_System__
#define __Poisson_CG_System__

#include <nova/Tools/Krylov_Solvers/Krylov_System_Base.h>
#include "../CG_Vector.h"
#include "../Convergence_Norm_Helper.h"
#include "Grid_Hierarchy_Projection.h"
#include "../Inner_Product_Helper.h"
#include "Combination_Multiply_Helper.h"

// #include "Poisson_Multigrid_Solver.h"

namespace Nova{
template<class Base_struct_type,class T,int d>
class Poisson_CG_System: public Krylov_System_Base<T>
{
    using Base                      = Krylov_System_Base<T>;
    using Vector_Base               = Krylov_Vector_Base<T>;
    using Channel_Vector            = Vector<T Base_struct_type::*,d>;
    using Hierarchy                 = Grid_Hierarchy<Base_struct_type,T,d>;
    using Hierarchy_Projection      = Grid_Hierarchy_Projection<Base_struct_type,T,d>;

  public:
    Hierarchy& hierarchy;
    T Base_struct_type::* alpha1_channel;
    T rho1,rho2;
    Poisson_CG_System(Hierarchy& hierarchy_input,T Base_struct_type::* alpha1_channel_input,T rho1_input,T rho2_input)
        :Base(false,false),hierarchy(hierarchy_input), alpha1_channel(alpha1_channel_input),rho1(rho1_input),rho2(rho2_input)
    {
    }

    ~Poisson_CG_System() {}

    void Set_Boundary_Conditions(Vector_Base& v) const {}
    void Project_Nullspace(Vector_Base& v) const {}

    void Multiply(const Vector_Base& v,Vector_Base& result) const
    {
        const int level=0;
        T Base_struct_type::* v_channel         = CG_Vector<Base_struct_type,T,d>::Cg_Vector(v).channel;
        T Base_struct_type::* result_channel    = CG_Vector<Base_struct_type,T,d>::Cg_Vector(result).channel;
        Combination_Multiply_Helper<Base_struct_type,T,d>(hierarchy,hierarchy.Allocator(level),hierarchy.Blocks(level),
                            v_channel,alpha1_channel,result_channel,rho1,rho2);
        
    }

    void Project(Vector_Base& v) const
    {
        T Base_struct_type::* v_channel         = CG_Vector<Base_struct_type,T,d>::Cg_Vector(v).channel;

        for(int level=0;level<hierarchy.Levels();++level)
            Clear_Non_Active<Base_struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v_channel);
    }

    double Inner_Product(const Vector_Base& v1,const Vector_Base& v2) const
    {
        const Hierarchy& v1_hierarchy           = CG_Vector<Base_struct_type,T,d>::Hierarchy(v1);
        const Hierarchy& v2_hierarchy           = CG_Vector<Base_struct_type,T,d>::Hierarchy(v2);
        T Base_struct_type::* const v1_channel  = CG_Vector<Base_struct_type,T,d>::Cg_Vector(v1).channel;
        T Base_struct_type::* const v2_channel  = CG_Vector<Base_struct_type,T,d>::Cg_Vector(v2).channel;
        assert(&hierarchy == &v1_hierarchy);
        assert(&hierarchy == &v2_hierarchy);

        double result=0;

        for(int level=0;level<hierarchy.Levels();++level)
            Inner_Product_Helper<Base_struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),v1_channel,
                                                  v2_channel,result,(unsigned)Cell_Type_Interior);

        return result;
    }

    T Convergence_Norm(const Vector_Base& v) const
    {
        T Base_struct_type::* v_channel         = CG_Vector<Base_struct_type,T,d>::Cg_Vector(v).channel;
        T max_value=(T)0.;

        for(int level=0;level<hierarchy.Levels();++level)
            Convergence_Norm_Helper<Base_struct_type,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                                                     v_channel,max_value,(unsigned)Cell_Type_Interior);

        return max_value;
    }

    void Apply_Preconditioner(const Vector_Base& r,Vector_Base& z) const
    {
    }
};
}
#endif
