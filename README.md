## Multicellular in-silico model reconstruction

# Description
This project presents a functional multicellular in-silico model reconstruction method. <br /> 
Its files were developed using a bespoke human genome constraint-based reconstruction workflow that <br />
combines aspects from an updated mCADRE & Metabotools algorithms, the novel redHuman algorithm, along with <br />
13C-metabolic flux analysis. 

The project aims to reproduce vital metabolic behaviours of in-vitro models in in-silico models, <br />
including cell growth predictions, respiration rates, and observations that suggest cross-shuttling <br />
of redox-active metabolites between cells.<br />

This repository contains: <br /> 
The mCADRE algorithm updated and optimised by Elias Vera-Siguenza. <br />
Originally Created by Wang et al. <br />
Wang, Y., Eddy, J.A. & Price, N.D. <br />
Reconstruction of genome-scale metabolic models for 126 human tissues using mCADRE. <br />
BMC Syst Biol 6, 153 (2012). https://doi.org/10.1186/1752-0509-6-1 <br />

The Multi-Objective Flux Analysis Library. <br />
Copyright (c) 2021, <br />
Lawrence Livermore National Security, LLC. <br />
Produced at the Lawrence Livermore National Laboratory. <br />
Written by Marc Griesemer and Ali Navid.

A compiled version of the Very Fast Flux Variability Analysis (VFFVA). <br />
File: veryfastFVA.c Version 0.4 <br />
Licence CC BY 4.0 : Free to share and modify. <br />
Author : Marouen Ben Guebila - marouen.benguebila@uni.lu <br />
based on Open MP/MPI parallel optimization of IBM's CPlex fastFVA. <br />

# How to Use
The in-silico model reconstruction files can be accessed and used by running the provided scripts <br />
along with the latest version of the human metabolic network matrix (Recon3D). <br />

# Requirements
MATLAB R2021a Update 3 (9.10.9.1684407) 64-bit (maci64) <br />
COnstraint-Based Reconstruction (COBRA - 2.39.2) <br />
IBM CPLEX_Studio <br />
MATLAB Bioinformatics Toolbox <br />
Affymetrix® CDF Data File Formats: The CDF file describes the layout for an Affymetrix GeneChip array. <br /> 
Parallel Computing Toolbox <br />
Obtain the Recon3D COBRA model MATLAB structure <br />
from: https://www.vmh.life/#downloadview <br />
Very Fast Variability Analysis library <br />
Note: <br />
IBM CPLEX has stopped its support for MATLAB since version 12.10, <br />
therefore VFFVA can be an alternative to access new CPLEX versions <br />
through its C API - instructions on how to achieve this can be found <br />
here: https://github.com/marouenbg/VFFVA <br />

The updated mCADRE script in this project was succesfully ran in less than 24 hours using a <br />
MacBook Pro 2.8 GHz Quad-Core Intel Core i7, <br />
Radeon Pro 555 2GB Intel HD Graphics 630 1536 MB GPU, <br />
16 GB 2133 MHz LPDDR3 RAM. (OSX: Ventura) <br />

Summary of system state at the time of solving:

			Support 	   LP 	 MILP 	   QP 	 MIQP 	  NLP 	   EP 
	------------------------------------------------------------------------------
	gurobi       	active        	    1 	    1 	    1 	    1 	    - 	    -
	ibm_cplex    	active        	    1 	    1 	    1 	    1 	    - 	    -
	tomlab_cplex 	active        	    0 	    0 	    0 	    0 	    - 	    -
	glpk         	active        	    1 	    1 	    - 	    - 	    - 	    -
	mosek        	active        	    0 	    - 	    0 	    - 	    - 	    0
	matlab       	active        	    1 	    - 	    - 	    - 	    1 	    -
	pdco         	active        	    1 	    - 	    1 	    - 	    - 	    1
	quadMinos    	active        	    1 	    - 	    - 	    - 	    - 	    -
	dqqMinos     	active        	    1 	    - 	    1 	    - 	    - 	    -
	cplex_direct 	active        	    0 	    0 	    0 	    - 	    - 	    -
	cplexlp      	active        	    1 	    - 	    - 	    - 	    - 	    -
	qpng         	passive       	    - 	    - 	    1 	    - 	    - 	    -
	tomlab_snopt 	passive       	    - 	    - 	    - 	    - 	    0 	    -
	lp_solve     	legacy        	    1 	    - 	    - 	    - 	    - 	    -
	------------------------------------------------------------------------------
	Total        	-             	    9 	    3 	    5 	    2 	    1 	    1

 + Legend: - = not applicable, 0 = solver not compatible or not installed, 1 = solver installed.<br />

# Contributors
Elias Vera-Siguenza (e.vera-siguenza@bham.ac.uk)

# MIT License:  
This document is based on the original mCADRE workflow. <br />
Briefly, the mCADRE is an algorithm aimed at reconstructing genome-scale  <br />
metabolic models for human tissues. It involves a data-driven approach  <br />
to create tissue-specific metabolic models. The algorithm utilizes gene  <br />
expression data, protein-protein interaction networks, and metabolic  <br />
pathway information to create tissue-specific metabolic models.  <br />
The resulting models can be used to predict metabolic fluxes,  <br />
simulate metabolic perturbations, and identify potential drug targets  <br />
specific to each tissue. The mCADRE approach allows for a comprehensive  <br />
understanding of tissue-specific metabolism. From its inception, it was  <br />
originally designed for use with the global knowledgebase of metabolic  <br />
functions categorized for the human genome (Human Recon 1).  <br />
Wang, Y., Eddy, J.A. & Price, N.D. Reconstruction of genome-scale  <br />
metabolic models for 126 human tissues using mCADRE.  <br />
BMC Syst Biol 6, 153 (2012). https://doi.org/10.1186/1752-0509-6-153 <br />

The purpose of this MATLAB library is to provide an updated version of  <br />
the original mCADRE. Our library has been developed to improve its  <br />
applicability and it has been optimised for its implementation with the  <br />
Recon3D metabolic network model (the latest iteration of Recon models).  <br />
This updated version also implements very fast flux variability analysis,  <br />
which allows for more efficient prediction of metabolic fluxes and  <br />
identification of metabolic targets. By incorporating the Recon3D model,  <br />
the updated mCADRE algorithm includes additional metabolic reactions and  <br />
pathways, leading to a more comprehensive analysis of metabolic networks  <br />
through a computationally less expensive algorithm; i.e., in less time.  <br />
Overall, the updated mCADRE algorithm is expected to optimise an already   <br />
accurate and efficient tool for predicting metabolic behavior and  <br />
identifying potential drug targets. <br />

The authors of this updated version aknowledge the original authors' work <br />
by merely stating that their contributions are just an update and not a <br />
means to change the original and itended outcome of the algorithm.  <br />
We thus adhere to the MIT License:  <br />

Copyright (c) 2023 Elias Vera-Siguenza <br />

Permission is hereby granted, free of charge, to any person obtaining a copy <br />
of the updated mCADRE version and associated documentation files  <br />
(the "Software"), to deal in the Software without restriction, including  <br />
without limitation the rights to use, copy, modify, merge, publish,  <br />
distribute, sublicense, and/or sell copies of the Software, and to permit  <br />
persons to whom the Software is furnished to do so, subject to the  <br />
following conditions: <br />

The above copyright notice and this permission notice shall be included  <br />
in all copies or substantial portions of the Software. <br />

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS  <br />
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, <br />
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL  <br />
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER <br />
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING  <br />
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER  <br />
DEALINGS IN THE SOFTWARE. <br />

By using, modifying, or distributing the updated mCADRE version,  <br />
Licensee acknowledges and agrees that Licensor shall not be liable for  <br />
any damages or losses, including but not limited to direct, indirect,  <br />
incidental, or consequential damages, arising from or related to the use  <br />
or inability to use the updated mCADRE version, even if Licensor has been  <br />
advised of the possibility of such damages. <br />

This license agreement does not create any form of partnership, agency,  <br />
or employment relationship between Licensor and Licensee. The license  <br />
granted herein is not transferable or assignable by Licensee without the  <br />
prior written consent of Licensor. <br />

By using, modifying, or distributing the updated mCADRE version,  <br />
Licensee acknowledges and agrees to indemnify and hold Licensor harmless  <br />
from any claims, damages, losses, or expenses, including reasonable  <br />
attorneys' fees, arising from Licensee's use or distribution of the  <br />
updated mCADRE version. <br />

This Agreement shall be governed by and construed in accordance with the  <br />
laws of United Kingdom, without regard to its conflict of law provisions. <br />
