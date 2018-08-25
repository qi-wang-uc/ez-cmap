#include <fstream>
#include <sstream>
#include "config.hpp"
#include "util.hpp"

bool Config::readConfig(const Str& inp_name) {
    std::cout << "ReadConfig> Reading configuration from file: " 
              << inp_name << std::endl;
    std::ifstream inp_file(inp_name);
    if(!inp_file.is_open()) 
        return error_exit<bool>("Cannot open input file.", false);
    Str each_line;
    Str buffer;
    StrQue commands;
    std::stringstream each_stream;
    
    while (getline(inp_file, each_line)) {
        if(each_line[0]=='#' || each_line.empty()) continue;
        each_stream.clear();
        each_stream.str(each_line);
        while (each_stream >> buffer) {
            commands.push(buffer);
        }
    }
    inp_file.close();
    while(!commands.empty()){
        if ("psfname"==commands.front()) {
            commands.pop();
            this->psf_name = commands.front();
            commands.pop();
        } else if ("dcdname"==commands.front()) {
            commands.pop();
            this->dcd_name = commands.front();
            commands.pop();
        } else if ("dcdspec"==commands.front()) {
            commands.pop();
            this->dcd_spec = commands.front();
            commands.pop();
        } else if ("cutoff" == commands.front()) {
            commands.pop();
            this->r_cutoff = std::stof(commands.front());
            commands.pop();
        } else if ("sele1"==commands.front()) {
            commands.pop();
            StrQue sele_str;
            while(commands.front()!="end") {
                sele_str.push(commands.front());
                commands.pop();
            }
            this->dmat_sele.first = retrieveSeleContext(sele_str);
        } else if ("sele2"==commands.front()) {
            commands.pop();
            StrQue sele_str;
            while(commands.front()!="end") {
                sele_str.push(commands.front());
                commands.pop();
            }
            this->dmat_sele.second = retrieveSeleContext(sele_str);
        } else if ("outpref"==commands.front()) {
            commands.pop();
            this->out_pref = commands.front();
        } else {
            commands.pop();
        }
    }
    return true;
}

Str Config::getConfigPsfName() const {
    return this->psf_name;
}

Str Config::getConfigDcdName() const {
    return this->dcd_name;
}

Str Config::getConfigDcdSpec() const {
    return this->dcd_spec;
}

float Config::getConfigCutoff() const {
    return this->r_cutoff;
}

std::pair<SeleContext, SeleContext> Config::getConfigSeleContext() const {
    return this->dmat_sele;
}

Str Config::getConfigOutPref() const {
    return this->out_pref;
}

void Config::printConfig() const {
    std::cout << "ReadConfig> After reading, the following configure will be used:" 
              << std::endl
              << "PSF NAME: " << this->psf_name << std::endl
              << "DCD NAME: " << this->dcd_name << std::endl
              << "R_CUTOFF: " << this->r_cutoff << std::endl
              << "DCD FRAMES TO USE: " << this->dcd_spec << std::endl
              << "SELECTION 1: " << std::endl
              << "  residue " << this->dmat_sele.first.resi_range.first 
              << " to "       << this->dmat_sele.first.resi_range.second 
              << std::endl
              << "  type " << printStrVec(this->dmat_sele.first.resi_type)
              << std::endl << "  segid " << this->dmat_sele.first.resi_segid 
              << std::endl
              << "SELECTION 2: " << std::endl
              << "  residue " << this->dmat_sele.second.resi_range.first 
              << " to "  << this->dmat_sele.second.resi_range.second
              << std::endl
              << "  type " << printStrVec(this->dmat_sele.second.resi_type)
              << std::endl << "  segid " << this->dmat_sele.second.resi_segid 
              << std::endl
              << "OUT PREF: " << this->out_pref 
    << std::endl;
}
