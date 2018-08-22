#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <regex>
#include <set>
#include <map>
#include <iterator>
#include <iomanip>
#include "ez-cmap.hpp"

std::vector<PsfAtom> PsfData;
std::vector<float> Dmatrix;
std::map<std::pair<size_t, size_t>, size_t> Cmap;
std::vector<float> dcd_buffer;

std::vector<size_t> AtomSele1;
std::vector<size_t> AtomSele2;
std::vector<Vec3d>  CordSele1;
std::vector<Vec3d>  CordSele2;
std::map<size_t, size_t> ResMapSele1;
std::map<size_t, size_t> ResMapSele2;

Config readConfig(const std::string& inp_name);
Range getRange(const std::string& inp_str);
Sele  getSele(std::queue<std::string> inp_str);
bool readPsfData(const std::string& inp_name);
bool procDcdData(const std::string& inp_name, const float& r_cutoff, const std::string& out_pref);
size_t getAtomCand(const Sele& sele, std::vector<size_t>& atomsele, std::map<size_t, size_t>& mapsele);
void initCoor(const size_t& dim1, const size_t& dim2);
void initCmap(void);

int main(int argc, char* argv[]) {
    if (1==argc) {
        std::cout << "ERROR> Missing input file." << std::endl;
        return 0;
    }
    std::string inp_name = argv[1];
    auto config = readConfig(inp_name);
    config.printConfig();
    if (!readPsfData(config.psf_name))
        return 0;
    const auto dim1 = getAtomCand(config.dmat_sele.first, AtomSele1, ResMapSele1);
    const auto dim2 = getAtomCand(config.dmat_sele.second, AtomSele2, ResMapSele2);
    initCoor(dim1, dim2);
    initCmap();
    if (!procDcdData(config.dcd_name, config.r_cutoff, config.out_pref))
        return 0;
    
    return 0;
}

Config readConfig(const std::string& inp_name) {
    Config result;
    std::ifstream inp_file(inp_name);
    std::string each_line;
    std::string buffer;
    std::stringstream each_stream;
    std::queue<std::string> commands;
    
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
            result.psf_name = commands.front();
            commands.pop();
        } else if ("dcdname"==commands.front()) {
            commands.pop();
            result.dcd_name = commands.front();
            commands.pop();
        } else if ("cutoff" == commands.front()) {
            commands.pop();
            result.r_cutoff = std::stof(commands.front());
            commands.pop();
        } else if ("sele1"==commands.front()) {
            commands.pop();
            std::queue<std::string> sele_str;
            while(commands.front()!="end") {
                sele_str.push(commands.front());
                commands.pop();
            }
            result.dmat_sele.first = getSele(sele_str);
        } else if ("sele2"==commands.front()) {
            commands.pop();
            std::queue<std::string> sele_str;
            while(commands.front()!="end") {
                sele_str.push(commands.front());
                commands.pop();
            }
            result.dmat_sele.second = getSele(sele_str);
        } else if ("outname"==commands.front()) {
            commands.pop();
            result.out_pref = commands.front();
        } else {
            commands.pop();
        }
    }
    return result;
}

Range getRange(std::string& inp_str) {
    Range result;
    auto sep_pos = inp_str.find(':');
    result.first = std::atoi(inp_str.substr(0, sep_pos).c_str());
    result.second = std::atoi(inp_str.substr(sep_pos+1).c_str());
    return result;
}

Sele getSele(std::queue<std::string> inp_str) {
    Sele result;
    while(!inp_str.empty()) {
        if ("resid"==inp_str.front()) {
            inp_str.pop();
            result.resi_range = getRange(inp_str.front());
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

bool readPsfData(const std::string& inp_name) {
    std::cout << "ReadPSF> Reading PSF info from file: " << inp_name << std::endl;
    std::ifstream inp_file(inp_name);
    if(!inp_file.is_open()) return false;
    std::regex r_psfatom("^\\s+\\d+\\s\\w+\\s\\d+\\s+\\w+\\s{2}\\w+.*");
    std::string each_line;
    std::stringstream each_stream;
    while (getline(inp_file, each_line)) {
        if (std::regex_match(each_line, r_psfatom)) {
            PsfAtom tmp_atom;
            each_stream.clear();
            each_stream.str(each_line);
            each_stream >> tmp_atom.atomid >> tmp_atom.segid >> tmp_atom.resid 
                        >> tmp_atom.resname >> tmp_atom.type;
            PsfData.push_back(tmp_atom);
        }
    }
    std::cout << "ReadPSF> After reading, (" << PsfData.size()
              << ") atoms were recorded." << std::endl;
    return true;
}

size_t getAtomCand(const Sele& sele, std::vector<size_t>& atomsele, std::map<size_t, size_t>& mapsele) {
    atomsele.resize(0);
    std::vector<PsfAtom> check_resid;
    std::vector<PsfAtom> check_type;
    std::vector<PsfAtom> check_segid;
    for (const auto& each_atom : PsfData) {
        if (std::find(sele.resi_type.begin(), sele.resi_type.end(), each_atom.type)!=sele.resi_type.end()) {
            check_type.push_back(each_atom);
        }    
        if (each_atom.segid==sele.resi_segid) {
            check_segid.push_back(each_atom);
        }
    }
    std::set_intersection(check_type.begin(), check_type.end(), 
        check_segid.begin(), check_segid.end(), std::back_inserter(check_resid));
    for (const auto& each_atom : check_resid) {
        if (each_atom.resid >= sele.resi_range.first && each_atom.resid <= sele.resi_range.second) {
            atomsele.push_back(each_atom.atomid);
            mapsele[each_atom.atomid] = each_atom.resid;
        }
    }
    return atomsele.size();
}

void initCoor(const size_t& dim1, const size_t& dim2) {
    std::cout << "Init> Initializing buffer for DMAT-related calculations" << std::endl;
    CordSele1.resize(dim1);
    CordSele2.resize(dim2);
    Dmatrix.resize(dim1*dim2);
    std::fill(Dmatrix.begin(), Dmatrix.end(), 0.0);
}

void initCmap(void) {
    std::cout << "Init> Initializing contact map" << std::endl;
    for(const auto& atom_i : AtomSele1) {
        auto resi_i = ResMapSele1[atom_i];
        for(const auto& atom_j : AtomSele2) {
            auto resi_j = ResMapSele2[atom_j];
            // std::cout << resi_i << " " << resi_j << std::endl;
            Cmap[std::make_pair(resi_i, resi_j)] = 0;
        }
    }
}

bool procDcdData(const std::string& inp_name, const float& r_cutoff, const std::string& out_pref) {
    std::cout << "ReadDCD> Processing DCD info from: " << inp_name << std::endl;
    std::ifstream inp_file(inp_name, std::ios::binary);
    if(!inp_file.is_open()) {
        std::cout << "ERROR> Cannot open DCD file." << std::endl;
        return false;
    }
    // Process DCD header
    DCD_Info dcd_info;
    char dcd_head1[100];
    char dcd_head2[16];
    char title[80];
    inp_file.read(dcd_head1, 100);
    dcd_info.dcd_header1.assign(dcd_head1, 100);
    if(strncmp(&dcd_head1[4], "CORD", 4)!=0) {
        std::cout << "ERROR> Wrong DCD format" << std::endl;
        return false;
    }
    int n_file = *(int*)(&dcd_head1[8]);
    dcd_info.n_frame = n_file;
	int n_priv = *(int*)(&dcd_head1[12]);
	int n_savc = *(int*)(&dcd_head1[16]);
	int n_step = *(int*)(&dcd_head1[20]);
	float delta = *(float*)(&dcd_head1[44]);
	int q_cell = *(int*)(&dcd_head1[48]);
    dcd_info.q_cell = q_cell;
	int c24tag = *(int*)(&dcd_head1[84]);
	if(c24tag!=24) std::cout << "ReadDCD> NOT NAMD trajectory format" << std::endl;
    std::cout << "ReadDCD> After reading, the following information were found:" << std::endl
              << "ReadDCD> NFILE=" << n_file << " NPRIV=" << n_priv
              <<" NSAVC=" << n_savc << " NSTEP=" << n_step << std::endl
              << "ReadDCD> DELTA=" << delta
              << (q_cell==1?"  PBC cells detected.":"  PBC cells NOT detected.") << std::endl;
    inp_file.read(title, 80);
    dcd_info.dcd_remark1.assign(title, 80);
	std::cout << "ReadDCD> " << title << std::endl;
	inp_file.read(title, 80);
    dcd_info.dcd_remark2.assign(title, 80);
	std::cout << "ReadDCD> " << title << std::endl;
	inp_file.read(dcd_head2, 16);
	int n_atom = *(int*)(&dcd_head2[8]);
    dcd_info.n_atom = n_atom;
	std::cout << "ReadDCD> (" << n_atom << ") atoms found in trajectory file." << std::endl;
    dcd_info.x_offset = (q_cell==1) ? 14          : 0;	
	dcd_info.y_offset = (q_cell==1) ? n_atom+16   : n_atom+2;
	dcd_info.z_offset = (q_cell==1) ? 2*n_atom+18 : 2*n_atom+4;
	dcd_info.sz_frame = (q_cell==1) ? (3*(4*n_atom+8)+56) : (3*(4*n_atom+8));
    dcd_buffer.resize(dcd_info.sz_frame);
    // Processing each frame, organized as follows: 
	// cell_offset + coor_pad + xcoor + 2*coor_pad + y_coor + 2*coor_pad + z_coor + coor_pad.
	// (14 or 0)   +    (1)      + n_atom+      (2)      + n_atom+      (2)      + natom  +    (1).
    const auto dim1 = CordSele1.size();
    const auto dim2 = CordSele2.size();
    const float rcut2 = r_cutoff*r_cutoff;
    for (auto iframe=0; iframe < dcd_info.n_frame; ++iframe) {
        std::cout << "ReadDCD> Processing frame: " << (iframe+1) << std::endl;
        inp_file.read(reinterpret_cast<char*>(&dcd_buffer[0]), dcd_info.sz_frame);
        for(size_t iatom=0; iatom<CordSele1.size(); ++iatom) {
		    CordSele1[iatom].x = dcd_buffer[dcd_info.x_offset + AtomSele1[iatom]];
            CordSele1[iatom].y = dcd_buffer[dcd_info.y_offset + AtomSele1[iatom]];
            CordSele1[iatom].z = dcd_buffer[dcd_info.z_offset + AtomSele1[iatom]];
	    }
        for(size_t iatom=0; iatom<CordSele2.size(); ++iatom) {
		    CordSele2[iatom].x = dcd_buffer[dcd_info.x_offset + AtomSele2[iatom]];
            CordSele2[iatom].y = dcd_buffer[dcd_info.y_offset + AtomSele2[iatom]];
            CordSele2[iatom].z = dcd_buffer[dcd_info.z_offset + AtomSele2[iatom]];
	    }
    /* Zero out Cmap*/
        for (auto& c : Cmap) {
            c.second = 0;
        }
    /* C++ code to build binary distance matrix */
        for(auto idx_i=0; idx_i<dim1; ++idx_i) {
            auto cord_i = CordSele1[idx_i];
            auto atom_i = AtomSele1[idx_i];
            auto resi_i = ResMapSele1[atom_i];
            for(auto idx_j=0; idx_j<dim2; ++idx_j) {
                auto cord_j = CordSele2[idx_j];
                auto atom_j = AtomSele2[idx_j];
                auto resi_j = ResMapSele2[atom_j];
                // Dmatrix[idx_i*dim1+idx_j] = atom_i.distsq(atom_j) < rcut2 ? 1 : 0;
                auto br_ij = cord_i.distsq(cord_j) < rcut2 ? 1 : 0;
                Cmap[std::make_pair(resi_i, resi_j)] += br_ij;
            }
        }
    /* Write to output */
        std::string out_name = out_pref + "-frame" + std::to_string(iframe+1) + ".dat";
        std::ofstream out_file(out_name);
        for (const auto& c : Cmap) {
            out_file << std::setw(5) << c.first.first << std::setw(5) << c.first.second << std::setw(5) << c.second << std::endl;
        }
        out_file.close();
    }
    
    return true;
}
