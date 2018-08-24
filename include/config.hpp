#ifndef CONFIG_HPP
#define CONFIG_HPP

#include "define.hpp"

class Config {
    protected:
        Str psf_name;
        Str dcd_name;
        float r_cutoff;
        std::pair<SeleContext, SeleContext> dmat_sele;
        Str out_pref;
    public:
    // setter
        bool readConfig(const Str& inp_name);
    // getter
        Str getConfigPsfName() const;
        Str getConfigDcdName() const;
        float getConfigCutoff() const;
        std::pair<SeleContext, SeleContext> getConfigSeleContext() const;
        Str getConfigOutPref() const;
        void printConfig() const;

};


#endif