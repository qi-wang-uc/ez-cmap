#ifndef PSF_HPP
#define PSF_HPP

#include <map>
#include "define.hpp"

/* Currently only [atom index], [segment name], [residue index], 
    [residue name] and [atom type] are supported */
struct PsfAtom {
    Int atomid;
    Str segid;
    Int resid;
    Str resname;
    Str type;
    PsfAtom() {}
    PsfAtom(Int atomid, Str segid, Int resid, Str resname, Str type):
        atomid(atomid),segid(segid),resid(resid),resname(resname),type(type) {}
    bool operator < (const PsfAtom& rhs) const {
        return atomid < rhs.atomid;
    }
};

class PSF {
    private:
        Int psf_natom;
    protected:
        std::vector<PsfAtom> PsfData;
    public:
        bool readPsfData(const Str& inp_name);
        std::vector<PsfAtom> const& getPsfData() const;
        Int  getPsfNatom(void) const;
};

#endif