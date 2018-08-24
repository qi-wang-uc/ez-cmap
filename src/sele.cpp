#include <iostream>
#include <algorithm>
#include "sele.hpp"
#include "util.hpp"

Selection::Selection(const SeleContext& sele_context, const std::vector<PsfAtom>& psf_data) {
    std::cout << "Selection> Creating selection .." << std::endl; 
    std::vector<PsfAtom> tmp_type;
    std::vector<PsfAtom> tmp_segid;
    std::vector<PsfAtom> tmp_resid;
    for (const auto& each_atom : psf_data) {
        if (std::find(sele_context.resi_type.begin(), sele_context.resi_type.end(), each_atom.type)!=sele_context.resi_type.end()) {
            tmp_type.push_back(each_atom);
        }    
        if (each_atom.segid==sele_context.resi_segid) {
            tmp_segid.push_back(each_atom);
        }
    }
    std::set_intersection(tmp_type.begin(), tmp_type.end(), 
        tmp_segid.begin(), tmp_segid.end(), std::back_inserter(tmp_resid));
    for (const auto& each_atom : tmp_resid) {
        if (each_atom.resid >= sele_context.resi_range.first && each_atom.resid <= sele_context.resi_range.second) {
            this->AtomIndex.push_back(each_atom.atomid);
            this->AtomToResMap[each_atom.atomid] = each_atom.resid;
        }
    }
    std::cout << "Selection> (" << AtomIndex.size() << ") atoms selected." << std::endl;
    this->sele_natom = AtomIndex.size();
}

Int Selection::getSeleNatom() const {
    return this->sele_natom;
}

std::vector<Int> const& Selection::getAtomIdVec() const {
    return this->AtomIndex;
}
Int Selection::getMappedRes(const Int& query) const {
    return this->AtomToResMap.at(query);
}
