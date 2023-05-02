## Multicellular in-silico model reconstruction

# Description
This project presents a functional multicellular in-silico model reconstruction method Its files were developed using a bespoke human genome constraint-based reconstruction workflow that combines aspects from an updated mCADRE & Metabotools algorithms, the novel redHuman algorithm, along with 13C-metabolic flux analysis. 

The project aims to reproduce vital metabolic behaviours of in-vitro models in in-silico models, including cell growth predictions, respiration rates, and observations that suggest cross-shuttling of redox-active metabolites between cells.

# How to Use
The in-silico model reconstruction files can be accessed and used by running the provided scripts aloong with the latest version of the human metabolic network matrix (Recon3D).

# Requirements
%% Requirements: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB R2021a Update 3 (9.10.9.1684407) 64-bit (maci64)
% COnstraint-Based Reconstruction (COBRA - 2.39.2) 
% IBM CPLEX_Studio
% Very Fast Variability Analysis library
% MATLAB Bioinformatics Toolbox
% Parallel Computing Toolbox
% Obtain the Recon3D COBRA model MATLAB structure 
% from: https://www.vmh.life/#downloadview

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: IBM CPLEX has stopped its support for MATLAB since version 12.10, 
% therefore VFFVA can be an alternative to access new CPLEX versions 
% through its C API - instructions on how to achieve this can be found
% herre: https://github.com/marouenbg/VFFVA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following script was succesfully ran in less than 24 hours using a 
% MacBook Pro 2.8 GHz Quad-Core Intel Core i7, 
% Radeon Pro 555 2GB Intel HD Graphics 630 1536 MB GPU, 
% 16 GB 2133 MHz LPDDR3 RAM. (OSX: Ventura)

%  Summary of available solvers and solver interfaces at the time of
%  running.
% 
% 			Support 	   LP 	 MILP 	   QP 	 MIQP 	  NLP 	   EP
% 	------------------------------------------------------------------------------
% 	gurobi       	active        	    1 	    1 	    1 	    1 	    - 	    -
% 	ibm_cplex    	active        	    1 	    1 	    1 	    1 	    - 	    -
% 	tomlab_cplex 	active        	    0 	    0 	    0 	    0 	    - 	    -
% 	glpk         	active        	    1 	    1 	    - 	    - 	    - 	    -
% 	mosek        	active        	    0 	    - 	    0 	    - 	    - 	    0
% 	matlab       	active        	    1 	    - 	    - 	    - 	    1 	    -
% 	pdco         	active        	    1 	    - 	    1 	    - 	    - 	    1
% 	quadMinos    	active        	    1 	    - 	    - 	    - 	    - 	    -
% 	dqqMinos     	active        	    1 	    - 	    1 	    - 	    - 	    -
% 	cplex_direct 	active        	    0 	    0 	    0 	    - 	    - 	    -
% 	cplexlp      	active        	    1 	    - 	    - 	    - 	    - 	    -
% 	qpng         	passive       	    - 	    - 	    1 	    - 	    - 	    -
% 	tomlab_snopt 	passive       	    - 	    - 	    - 	    - 	    0 	    -
% 	lp_solve     	legacy        	    1 	    - 	    - 	    - 	    - 	    -
% 	------------------------------------------------------------------------------
% 	Total        	-             	    9 	    3 	    5 	    2 	    1 	    1
% + Legend: - = not applicable, 0 = solver not compatible or not installed, 1 = solver installed.

# Contributors
Elias Vera-Siguenza (e.vera-siguenza@bham.ac.uk)

# License
This project is licensed under the MIT License - see the LICENSE.md file for details.
