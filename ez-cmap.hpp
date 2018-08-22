#ifndef EZDMAT_HPP
#define EZDMAT_HPP

#include <iostream>
#include <vector>
#include <utility>
#include <cmath>

using Range = std::pair<size_t, size_t>;

struct PsfAtom {
    size_t atomid;
    std::string segid;
    size_t resid;
    std::string resname;
    std::string type;
    PsfAtom() {}
    PsfAtom(size_t atomid, std::string segid, size_t resid, std::string resname, std::string type):
        atomid(atomid),segid(segid),resid(resid),resname(resname),type(type) {}
    bool operator < (const PsfAtom& rhs) const {
        return atomid < rhs.atomid;
    }
};

struct Vec3d {
    double x;
    double y;
    double z;
    double dist (const Vec3d& rhs) {
        return sqrt((x-rhs.x)*(x-rhs.x) + (y-rhs.y)*(y-rhs.y) + (z-rhs.z)*(z-rhs.z));
    }
    double distsq (const Vec3d& rhs) {
        return (x-rhs.x)*(x-rhs.x) + (y-rhs.y)*(y-rhs.y) + (z-rhs.z)*(z-rhs.z);
    }
};

struct Sele {
    Range resi_range;
    std::vector<std::string> resi_type;
    std::string resi_segid;
};

struct Config {
    std::string psf_name;
    std::string dcd_name;
    float r_cutoff;
    std::pair<Sele, Sele> dmat_sele;
    std::string out_pref;
    void printConfig(void);
};

void Config::printConfig(void) {
    std::cout << "After reading, the following configure will be used:" << std::endl;
    std::cout << "PSF NAME: " << this->psf_name << std::endl
              << "DCD NAME: " << this->dcd_name << std::endl
              << "R_CUTOFF: " << this->r_cutoff << std::endl;
    std::cout << "SELECTION 1: " << std::endl
              << "  residue " << this->dmat_sele.first.resi_range.first << " to " 
              << this->dmat_sele.first.resi_range.second << std::endl
              << "  type ";
              for(const auto& str : this->dmat_sele.first.resi_type)
                std::cout << str << " ";
    std::cout << std::endl << "  segid " << this->dmat_sele.first.resi_segid 
              << std::endl;
    std::cout << "SELECTION 2: " << std::endl
              << "  residue " << this->dmat_sele.second.resi_range.first << " to " 
              << this->dmat_sele.second.resi_range.second << std::endl
              << "  type ";
              for(const auto& str : this->dmat_sele.second.resi_type)
                std::cout << str << " ";
    std::cout << std::endl << "  segid " << this->dmat_sele.second.resi_segid 
              << std::endl;
    std::cout << "OUT PREF: " << this->out_pref << std::endl;
}

struct DCD_Info {
    // Header info, retrieved when reading header (except dcd_header2).
    std::string dcd_header1;    // [100] "CORD" and info of dcd stats.
    std::string dcd_remark1;    // [80] "REMARK CREATED BY ...
    std::string dcd_remark2;    // [80] "REMARK" + $DATE + $USER 
    std::string dcd_header2;    // [16] NATOM padded by dummy size_ts on both sides

    // Obtained from header.
    size_t n_atom  = 0;
    size_t n_frame = 0;
    size_t q_cell  = 0;
    
    // Caclulated based on above.
    size_t x_offset = 0;
	size_t y_offset = 0;
	size_t z_offset = 0;
	size_t sz_frame = 0;
};

#endif