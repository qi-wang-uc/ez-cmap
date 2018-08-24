#include <iostream>
#include "define.hpp"
#include "config.hpp"
#include "psf.hpp"
#include "sele.hpp"
#include "cmap.hpp"
#include "util.hpp"

int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);
    printProgName();
    if (1==argc) return error_exit<Int>("Missing input file", 1);

    // 1. Process input config
    Str config_name = argv[1];
    Config config;
    if(!config.readConfig(config_name)) return 1;
    config.printConfig();

    // 2. Read PSF and generate selection based on config
    PSF psf;
    if(!psf.readPsfData(config.getConfigPsfName())) return 1;

    // 3. Parse selection context to retrieve indices of selected atoms
    Selection sele1(config.getConfigSeleContext().first,  psf.getPsfData());
    Selection sele2(config.getConfigSeleContext().second, psf.getPsfData());

    // 4. Build and write contact map while processing DCD frame-wise.
    ContactMap cmap(sele1, sele2);
    if(!cmap.processDcdData(config, psf.getPsfNatom())) return 1;
    
    return 0;
}