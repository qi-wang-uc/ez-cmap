# ez-cmap
Calculate contact number map between selected residues

#### Usage:
- Compile with `make` in the program root folder
- Run with `./ez-cmap.exe input.dat`

#### Example Input:
```
# Read input
psfname input.psf
dcdname input.dcd
cutoff  8.0

# Define selections
sele1
    resid 1:340
    type CA N C O
    segid SEGI
end

sele2
    resid 1:340
    type CA N C O
    segid SEGI
end

# Outname prefix of cmap file, full name will be ${outname}-frame${framenum}.dat
outname cmap_out
```
*Note*: 
- The `resid` selection should be a complete string (e.g. `1:340`) without any space seperations.
- Currently the program only supports single residue range or segid but multiple atom types.

#### Example:
![example](test/test.png)

