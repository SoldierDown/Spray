//!#####################################################################
//! \file Poisson_Data.h
//!#####################################################################
// Class Poisson_Data
//######################################################################
#ifndef __Poisson_Data__
#define __Poisson_Data__

#include <stdint.h>

namespace Nova{
template<class T,class T_FLAGS=uint32_t>
struct Poisson_Data
{
    typedef T_FLAGS Flags_type;

    T_FLAGS flags;
    T ch0;      // X-face_velocity
    T ch1;      // Y-face_velocity
    T ch2;      // Z-face_velocity
    T ch3;      // X-Qc
    T ch4;      // Y-Qc
    T ch5;      // Z-Qc
    T ch6;      // density
    T ch7;      // density_backup
    T ch8;      // temp
    T ch9;      // temp
    T ch10;     // temp
    T ch11;     // temp
    T ch12;     // temp
    T ch13;     // temp
    T ch14;     // temp
    T ch15;     // temp
    T ch16;     // temp
    T ch17;     // temp
    T ch18;     // temp
    T ch19;     // temp
    T ch20;     // temp
    T ch21;     // temp
    T ch22;     // temp
    T ch23;     // temp
    T ch24;     // temp
};
}
#endif
