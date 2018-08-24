#ifndef SELE_HPP
#define SELE_HPP

#include <map>
#include "psf.hpp"
#include "define.hpp"

class Selection {
    private:
        Int sele_natom;
    protected:
        std::vector<Int> AtomIndex;
        std::map<Int, Int> AtomToResMap;
    public:
        Selection() {}
        Selection(const SeleContext& sele_context, const std::vector<PsfAtom>& psf_data);
        Int getSeleNatom() const;
        std::vector<Int> const& getAtomIdVec() const;
        Int getMappedRes(const Int& query) const;

};

#endif