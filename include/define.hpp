#ifndef DEFINE_HPP
#define DEFINE_HPP

#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <string>
#include <queue>

using Int   = size_t;
using Int2d = std::pair<size_t, size_t>;
using IntVec = std::vector<size_t>;

using Str = std::string;
using StrVec = std::vector<std::string>;
using StrQue = std::queue<std::string>;

struct Vec3d {
    float x;
    float y;
    float z;
    float dist (const Vec3d& rhs) {
        return sqrt((x-rhs.x)*(x-rhs.x) + (y-rhs.y)*(y-rhs.y) + (z-rhs.z)*(z-rhs.z));
    }
    float distsq (const Vec3d& rhs) {
        return (x-rhs.x)*(x-rhs.x) + (y-rhs.y)*(y-rhs.y) + (z-rhs.z)*(z-rhs.z);
    }
};

struct DCD_Info {
    // Header info, retrieved when reading header (except dcd_header2).
    Str dcd_header1;    // [100] "CORD" and info of dcd stats.
    Str dcd_remark1;    // [80] "REMARK CREATED BY ...
    Str dcd_remark2;    // [80] "REMARK" + $DATE + $USER 
    Str dcd_header2;    // [16] NATOM padded by dummy size_ts on both sides

    // Obtained from header.
    Int n_file = 0;
    Int n_priv  = 0;
    Int n_savc  = 0;
    Int n_step  = 0;
    Int n_atom  = 0;
    bool q_cell  = false;
    bool q_namd  = false;
    float delta = 0.0;
    
    // Caclulated based on above.
    Int x_offset = 0;
	Int y_offset = 0;
	Int z_offset = 0;
	Int sz_frame = 0;
};

struct SeleContext {
    Int2d resi_range;
    StrVec resi_type;
    Str resi_segid;
};

#endif