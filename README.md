## Multicellular in-silico model reconstruction

# Description
This project presents a functional multicellular in-silico model reconstruction method. <br /> 
Its files were developed using a bespoke human genome constraint-based reconstruction workflow that <br />
combines aspects from an updated mCADRE & Metabotools algorithms, the novel redHuman algorithm, along with <br />
13C-metabolic flux analysis. 

The project aims to reproduce vital metabolic behaviours of in-vitro models in in-silico models, <br />
including cell growth predictions, respiration rates, and observations that suggest cross-shuttling <br />
of redox-active metabolites between cells.

# How to Use
The in-silico model reconstruction files can be accessed and used by running the provided scripts <br />
along with the latest version of the human metabolic network matrix (Recon3D). <br />

# Requirements
MATLAB R2021a Update 3 (9.10.9.1684407) 64-bit (maci64) <br />
COnstraint-Based Reconstruction (COBRA - 2.39.2) <br />
IBM CPLEX_Studio <br />
MATLAB Bioinformatics Toolbox <br />
Parallel Computing Toolbox <br />
Obtain the Recon3D COBRA model MATLAB structure <br />
from: https://www.vmh.life/#downloadview <br />
Very Fast Variability Analysis library <br />
Note: <br />
IBM CPLEX has stopped its support for MATLAB since version 12.10, <br />
therefore VFFVA can be an alternative to access new CPLEX versions <br />
through its C API - instructions on how to achieve this can be found <br />
herre: https://github.com/marouenbg/VFFVA <br />

The following script was succesfully ran in less than 24 hours using a <br />
MacBook Pro 2.8 GHz Quad-Core Intel Core i7, <br />
Radeon Pro 555 2GB Intel HD Graphics 630 1536 MB GPU, <br />
16 GB 2133 MHz LPDDR3 RAM. (OSX: Ventura) <br />


# Contributors
Elias Vera-Siguenza (e.vera-siguenza@bham.ac.uk)

# License
This project is licensed under the MIT License - see the LICENSE.md file for details.
