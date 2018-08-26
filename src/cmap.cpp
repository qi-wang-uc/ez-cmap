#include <fstream>
#include <cstring>
#include <iomanip>
#include "cmap.hpp"
#include "define.hpp"
#include "util.hpp"

ContactMap::ContactMap(Selection sele1, Selection sele2) {
    std::cout << "CMAP> Initializing contact map ..." << std::endl;
    this->AtomSele1 = std::move(sele1);
    this->AtomSele2 = std::move(sele2);
    for(const auto& atom_i : this->AtomSele1.getAtomIdVec()) {
        auto resi_i = this->AtomSele1.getMappedRes(atom_i);
        for(const auto& atom_j : this->AtomSele2.getAtomIdVec()) {
            auto resi_j = this->AtomSele2.getMappedRes(atom_j);
            Cmap[std::make_pair(resi_i, resi_j)] = 0;
        }
    }
    std::cout << "CMAP> Initializing coordinates for selection ..." << std::endl;
    this->CordSele1.resize(sele1.getSeleNatom());
    this->CordSele2.resize(sele2.getSeleNatom());
}

bool ContactMap::processDcdData(const Config& config, const Int& psf_natom) {
    auto dcd_name = config.getConfigDcdName();
    auto out_pref = config.getConfigOutPref();
    auto r_cutoff = config.getConfigCutoff();
    auto dcd_spec = config.getConfigDcdSpec();
    // Sanity check
    std::cout << "ReadDCD> Processing DCD info from: " << dcd_name << std::endl;
    std::ifstream inp_file(dcd_name, std::ios::binary);
    if(!inp_file.is_open())
        return error_exit<bool>("Cannot open DCD file.", false);
    // Process DCD header
    DCD_Info dcd_info;
    char dcd_head1[100];
    char dcd_head2[16];
    char title[80];
    inp_file.read(dcd_head1, 100);
    dcd_info.dcd_header1.assign(dcd_head1, 100);
    if(strncmp(&dcd_head1[4], "CORD", 4)!=0)
        return error_exit<bool>("Wrong DCD format", false);
    dcd_info.n_file = *(int*)(&dcd_head1[8]);
    dcd_info.n_priv = *(int*)(&dcd_head1[12]);
    dcd_info.n_savc = *(int*)(&dcd_head1[16]);
    dcd_info.n_step = *(int*)(&dcd_head1[20]);
    dcd_info.delta = *(float*)(&dcd_head1[44]);
    dcd_info.q_cell = (1==*(int*)(&dcd_head1[48])) ? true:false;
    dcd_info.q_namd = (24==*(int*)(&dcd_head1[84])) ? true:false;
    inp_file.read(title, 80);
    dcd_info.dcd_remark1.assign(title, 80);
	inp_file.read(title, 80);
    dcd_info.dcd_remark2.assign(title, 80);
	inp_file.read(dcd_head2, 16);
    dcd_info.n_atom = *(int*)(&dcd_head2[8]);

    printDcdInfo(dcd_info);
    if(dcd_info.n_atom != psf_natom) {
        return error_exit<bool>("Mismatch of Natom in PSF and DCD.", false);
    }
    // Update offset information
    dcd_info.x_offset = (dcd_info.q_cell) ? 14 : 0;	
	dcd_info.y_offset = (dcd_info.q_cell) ? dcd_info.n_atom+16 : dcd_info.n_atom+2;
	dcd_info.z_offset = (dcd_info.q_cell) ? 2*dcd_info.n_atom+18 : 2*dcd_info.n_atom+4;
	dcd_info.sz_frame = (dcd_info.q_cell) ? (3*(4*dcd_info.n_atom+8)+56) : (3*(4*dcd_info.n_atom+8));
    // Process dcd frame request
    auto nframe = atoi(dcd_spec.c_str());
    bool is_single_frame = true;
    if(0==nframe && dcd_spec!="all" || nframe < 0) {
        return error_exit("Invalid dcd frame specification.", false);
    } else if (0==nframe && dcd_spec=="all") {
        is_single_frame = false;
        nframe = dcd_info.n_file;
    } else if (nframe > dcd_info.n_file) {
        std::cout << "ReadDCD> Request frame not available, last frame will be used." << std::endl;
        is_single_frame = true;
        nframe = dcd_info.n_file;
    }
    // Process each frame
    std::vector<float> dcd_buffer(dcd_info.sz_frame);
    const auto dim1 = this->CordSele1.size();
    const auto dim2 = this->CordSele2.size();
    const float rcut2 = r_cutoff*r_cutoff;
    for (auto iframe=0; iframe < nframe; ++iframe) {
        if(is_single_frame && (iframe+1)!=nframe) continue;
        std::cout << "ReadDCD> Processing frame: " << (iframe+1) << std::endl;
        inp_file.read(reinterpret_cast<char*>(&dcd_buffer[0]), dcd_info.sz_frame);
        for(size_t iatom=0; iatom<CordSele1.size(); ++iatom) {
		    CordSele1[iatom].x = dcd_buffer[dcd_info.x_offset + this->AtomSele1.getAtomIdVec()[iatom]];
            CordSele1[iatom].y = dcd_buffer[dcd_info.y_offset + this->AtomSele1.getAtomIdVec()[iatom]];
            CordSele1[iatom].z = dcd_buffer[dcd_info.z_offset + this->AtomSele1.getAtomIdVec()[iatom]];
	    }
        for(size_t iatom=0; iatom<CordSele2.size(); ++iatom) {
		    CordSele2[iatom].x = dcd_buffer[dcd_info.x_offset + this->AtomSele2.getAtomIdVec()[iatom]];
            CordSele2[iatom].y = dcd_buffer[dcd_info.y_offset + this->AtomSele2.getAtomIdVec()[iatom]];
            CordSele2[iatom].z = dcd_buffer[dcd_info.z_offset + this->AtomSele2.getAtomIdVec()[iatom]];
	    }
        buildContactMap(dim1, dim2, rcut2);
        writeContactMap(out_pref, iframe);
    }
    return true;
}

void ContactMap::buildContactMap(const Int& dim1, const Int& dim2, const float& rcut2) {
    for (auto& c : this->Cmap) {
            c.second = 0;
        }
    for(auto idx_i=0; idx_i<dim1; ++idx_i) {
        auto cord_i = this->CordSele1[idx_i];
        auto atom_i = this->AtomSele1.getAtomIdVec()[idx_i];
        auto resi_i = this->AtomSele1.getMappedRes(atom_i);
        for(auto idx_j=0; idx_j<dim2; ++idx_j) {
            auto cord_j = this->CordSele2[idx_j];
            auto atom_j = this->AtomSele2.getAtomIdVec()[idx_j];
            auto resi_j = this->AtomSele2.getMappedRes(atom_j);
            // Dmatrix[idx_i*dim1+idx_j] = atom_i.distsq(atom_j) < rcut2 ? 1 : 0;
            auto br_ij = cord_i.distsq(cord_j) < rcut2 ? 1 : 0;
            this->Cmap[std::make_pair(resi_i, resi_j)] += br_ij;
        }
    }
}

void ContactMap::writeContactMap(const Str& out_pref, const Int& iframe) {
    Str out_name = out_pref + "-frame" + std::to_string(iframe+1) + ".dat";
    std::ofstream out_file(out_name);
    for (const auto& c : this->Cmap) {
        out_file << std::setw(20) << c.first.first 
                 << std::setw(20) << c.first.second 
                 << std::setw(20) << c.second 
                 << std::endl;
    }
    out_file.close();
}