#include <iostream>
#include "util.hpp"

Str printStrVec(const StrVec& str_vec) {
    Str result;
    for(const auto& str : str_vec) {
        result += (str + " ");
    }
    return result;
}

Int2d retrieveRange(Str& inp_str) {
    Int2d result;
    Int sep_pos = inp_str.find(':');
    result.first = std::atoi(inp_str.substr(0, sep_pos).c_str());
    result.second = std::atoi(inp_str.substr(sep_pos+1).c_str());
    return result;
}

SeleContext retrieveSeleContext(StrQue inp_str) {
    SeleContext result;
    while(!inp_str.empty()) {
        if ("resid"==inp_str.front()) {
            inp_str.pop();
            result.resi_range = retrieveRange(inp_str.front());
        } else if ("type"==inp_str.front()) {
            inp_str.pop();
            while("segid"!=inp_str.front()) {
                result.resi_type.push_back(inp_str.front());
                inp_str.pop();
            }
        } else if ("segid"==inp_str.front()) {
            inp_str.pop();
            result.resi_segid = inp_str.front();
        } else {
            inp_str.pop();
        }
    }
    return result;
}

void printProgName() {
    std::cout << "EZCMAP> Generating residue contact number map. (version 1.1)"
              << std::endl;
}

void printDcdInfo(const DCD_Info& dcd_info) {
    std::cout << "ReadDCD> Summary of DCD information:" << std::endl
              << dcd_info.dcd_remark1 << std::endl
              << dcd_info.dcd_remark2 << std::endl
              << "ReadDCD> NFILE=" << dcd_info.n_file 
              << " NPRIV=" << dcd_info.n_priv
              << " NSAVC=" << dcd_info.n_savc 
              << " NSTEP=" << dcd_info.n_step << std::endl
              << "ReadDCD> DELTA=" << dcd_info.delta
              << "  PBC cells " << (dcd_info.q_cell ? "":"NOT ") << "detected" 
              << std::endl
              << "ReadDCD> DCD file is " << (dcd_info.q_namd ? "":"NOT ") << "NAMD format" 
              << std::endl
              << "ReadDCD> (" << dcd_info.n_atom 
              << ") atoms found in trajectory file." 
              << std::endl;
}