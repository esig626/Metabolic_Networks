%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Elias Siguenza
% Date: 02.05.2023
% The Tennant Lab
% Address: Institute of Metabolism and Systems Research, 
% College of Medical and Dental Sciences, University of Birmingham, 
% B15 2TT, England, United Kingdom.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ ME: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This document is based on the original mCADRE workflow.
% Briefly, the mCADRE is an algorithm aimed at reconstructing genome-scale 
% metabolic models for human tissues. It involves a data-driven approach 
% to create tissue-specific metabolic models. The algorithm utilizes gene 
% expression data, protein-protein interaction networks, and metabolic 
% pathway information to create tissue-specific metabolic models. 
% The resulting models can be used to predict metabolic fluxes, 
% simulate metabolic perturbations, and identify potential drug targets 
% specific to each tissue. The mCADRE approach allows for a comprehensive 
% understanding of tissue-specific metabolism. From its inception, it was 
% originally designed for use with the global knowledgebase of metabolic 
% functions categorized for the human genome (Human Recon 1). 
% Wang, Y., Eddy, J.A. & Price, N.D. Reconstruction of genome-scale 
% metabolic models for 126 human tissues using mCADRE. 
% BMC Syst Biol 6, 153 (2012). https://doi.org/10.1186/1752-0509-6-153
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% The purpose of this MATLAB library is to provide an updated version of 
% the original mCADRE. Our library has been developed to improve its 
% applicability and it has been optimised for its implementation with the 
% Recon3D metabolic network model (the latest iteration of Recon models). 
% This updated version also implements very fast flux variability analysis, 
% which allows for more efficient prediction of metabolic fluxes and 
% identification of metabolic targets. By incorporating the Recon3D model, 
% the updated mCADRE algorithm includes additional metabolic reactions and 
% pathways, leading to a more comprehensive analysis of metabolic networks 
% through a computationally less expensive algorithm; i.e., in less time. 
% Overall, the updated mCADRE algorithm is expected to optimise an already  
% accurate and efficient tool for predicting metabolic behavior and 
% identifying potential drug targets.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Requirements: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%% Initialise Workdesk

% Preload PATH and Initialise COBRA and IBM - CPLEX.
clear; close all; initCobraToolbox(false); changeCobraSolver('ibm_cplex');
filesFolder = '/Your/Path/To/mat/FILES';
mpsFolder = '/Your/Path/To/.MPS/Files';
CELFiles = '/Your/Path/To/.CEL/FILES';
clear ans; clc; cd(filesFolder)
load indGeneRxn_3 ; load rxn_map_matrix_3; load confidenceScores_3; 
load precursorMets_3; cd ..

%% Step 1: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read RNA data and set the Threshold defining what is usable and
% what is not. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For this section, the HG-U133_Plus_2.cdf file is needed. 
% Should be downloaded from Affymetrix website (Files provided)
% The next command will read the transcriptomic datasets I will use for the
% GEM reconstruction. They should be in .CEL format. All of which based 
% on Affymetrix Human Genome U133 Plus 2.0 Array platform.

% Open folder containing .CEL files:
cd(CELFiles)
% Run celintensityread to create RawStracture object. This command 
% generates a MATLAB structure field with several entries arranging the 
% transcriptomic datasets.
RawStructure = celintensityread('*', 'HG-U133_Plus_2.cdf','PMOnly',false);
clc
% Return to main folder: 
cd ..

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mas5Calls (Set the threshold for expressed and non-expressed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs the Wilcoxon signed rank-based gene expression presence/absence 
% detection algorithm first implemented in the Affymetrix Microarray Suite 
% V5. This specific algorithm uses a statistical test called the Wilcoxon 
% signed-rank test to help determine whether a gene is active (present) or 
% inactive (absent) in the sample. Usage: 
% [CALLS] = MAS5CALLS(RAW_STRUCT, TAU, ALPHA1, ALPHA2)
% RAW_STRUCT: 
% A struct of raw CEL data such as that returned by celintensityread
% TAU:		The tau cutoff for significant difference between PIi and MMi. 
% Defaults to 0.015.
% ALPHA1:	The p-value below which probesets are called 'present'.  
% Defaults to 0.04.
% ALPHA2:	The p-value above which probesets are called 'absent'.  
% Defaults to 0.06. 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% NOTE: 
% If parallel Toolbox is not available use deprecated mas5Calls 
% (original mCadre) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

calls = mas5Calls_new(RawStructure); 
dataset_probe=RawStructure.ProbeSetIDs;

%% Step 2:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mapping Genetic Expression data onto model reactions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load Affymetrix gene notation to Entrez gene notation 
% (manually curated to match those in Recon3D model - File Provided)
probes2idFile = fullfile(filesFolder, 'probes2id.mat');
load(probes2idFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the Recon3D model (curated with Entrez id's). This structure is 
% also equipped with a logic vector detailing reaction reversibility. 
% This can easily be obtained through indexing the model's lower bound 
% reaction values.
recon3DFile = fullfile(filesFolder, 'Recon3D.mat');
load(recon3DFile);
originalModel=model;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the gene-to-reaction associations included in Recon3D
% Optimised function. The original file intended for use with
%  Recon 1 had a numerical vector of Entrez gene id's. Recon3D has a cell
% array consisting of a string character of Entrez id's. I've taken this
% and converted it into a numeric vector to emulate the original procedure.
% This is a subcategory in the model object named genes2 in Recon3D.mat
% file. Using this file, now we are able to calculate the ubiquity scores 
% of metabolic genes at a considerably faster speed.
[parsedGPR,corrRxn] = extractGPRs_N(model);

parsedGPR = reshape(cellfun(@(x) x{1},...
                      parsedGPR, 'UniformOutput', false), size(parsedGPR));

[expression,expressionid]=mas5callToExpression(calls,...
    dataset_probe,probes,probe2geneid,model.genes2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We map the gene ubiquity scores into the model's metabolic networks to
% obtain expression-based evidence score for the Recon3D reactions:
[geneAssociatedRxns,reactionScores]= ...
                                geneExpression2rxnScore_N(expressionid,...
                                             expression,parsedGPR,corrRxn);


%% Step 3: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify core reactions. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We convert the numerical entrez gene ids to a cell array of strings
model.genes2=cellstr(num2str(model.genes2));

% Assign expression-based evidence to gene-associated reactions
numrxns=numel(model.rxns);
reactionScoresAll=zeros(numrxns,1);

% Add one one to reactionScores
reactionScoresAll(indGeneRxn_3,1)=reactionScores(:,1);
reactionScoresAll(isnan(reactionScoresAll))=0;

% Core reactions defined as reactions present in at least 50% of samples
% based on data. Can be a different threshold.
core=model.rxns(reactionScoresAll>=0.5);
 
% Check model consistency.
% Identify core reactions that can also carry flux
% If VFFA is not an option, you can run
% oldschool FVA function with a huge slowdown (Weeks). 
% You can alternatively skip and only check consistency of core reactions 
% using Jerby et al. method. (set fourth input to 0)

% To use Very Fast Flux Variability Analysis we set func input four : 1.
inactive_g=checkModelConsistency_N(model,[],[],1,mpsFolder);

% Find the Core Active Reactions:
active_g=setdiff(model.rxns,inactive_g);
core_active_g=intersect(core,active_g);

%% Step 4: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the connection-based evidence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind2active=ismember(model.rxns,active_g);
active_g=model.rxns(ind2active);
rxn_map_matrix_3=rxn_map_matrix_3-eye(size(rxn_map_matrix_3));
connect2core=rxn_map_matrix_3(ind2active,ind2active);
connect2core=connect2core./repmat(sum(connect2core,2),1,size(connect2core,2));
connect2coreScore=connect2core*reactionScoresAll(ind2active,1);
connect2coreScore(isnan(connect2coreScore))=0;

%% Step 5: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distinguish gene associated reaction NOT expressed (0 score)
% and non-gene associated reaction (0 score). In removal,gene
% associated reaction NOT expressed to be removed first. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

non_core=setdiff(model.rxns,union(core,setdiff(model.rxns,active_g)));

[c,ind]=ismember(non_core,active_g);

connect=connect2coreScore(ind(c));

[c,ind]=ismember(non_core,model.rxns);

exp_score=reactionScoresAll(ind(c),1);

confidence=confidenceScores_3(ind(c),1);

generxns=model.rxns(indGeneRxn_3);

ind2GeneNoExp=ismember(non_core,generxns);

exp_score(ind2GeneNoExp,1)=exp_score(ind2GeneNoExp,1)-1e-6;

%% Step 6: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort non_core rxns first by connectivity to core, then by expression, then
% by confidenceScores in Recon3D.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[c,indNonCore]=sortrows([connect,exp_score,confidence],[2,1,3]);

% Map connection-based,confidence score-based
% expression-based scores to non-core reactions and rank them

rxns2remove=non_core(indNonCore);

%% Step 7: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The initial tissue_model will be a model with all the 
% flux-carrying reactions in Recon3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu = 0;
count=0;

tissue_model=removeRxns(model,setdiff(model.rxns,active_g));
precursorMets=[precursorMets;nonEssentialAA;nucleotide;lipid];


%% Step 8:

if checkModelFunction(tissue_model,precursorMets)
disp('genericModel passed precursor metabolites test!');
%check if the generic model can pass the metabolite production test
core_removed=[];
%Record all the core reactions  removed during model pruning

while ~isempty(rxns2remove)
    
    rxn2remove=rxns2remove(1);
    disp('Next Reaction!');
    % Reactions that are inactivated by the removal of the current non-core reaction.
    inactive=checkModelConsistency_N(tissue_model,rxn2remove,core_active_g,1,mpsFolder);
    
    % Remove null-test reaction: (To avoid this, amend the VFFVA C++ code)
    inactive{1,end}=[];
    isEmpty = cellfun(@isempty, inactive);
    inactive(isEmpty) = [];

    disp('Done')
    
    if exp_score(ismember(non_core,rxn2remove))==-1e-6 &&... % If the non-core reaction is not present in ANY tissue samples
    numel(intersect(inactive,core_active_g))/numel(intersect(inactive,non_core))<=(1/3) &&... % If the number of non-core reactions being removed exceeds the number of core reactions being removed exceeds 3 to 1
    checkModelFunction(removeRxns(tissue_model,inactive),precursorMets)&&checkSalvage(removeRxns(tissue_model,inactive)) % And if removing these reactions won't affect tissue model's core functinalities
    
    % CheckSalvage tests the ability of the tissue model to use the 
    % salvage pathway to make purines, as de novo purine synthesis 
    % primarily occurs in the liver.  Not useful when the tissue is 
    % known to make purines using de novo synthesis pathway.
    
        tissue_model=removeRxns(tissue_model,inactive);
        rxns2remove=rxns2remove(~ismember(rxns2remove,inactive),1);
        display(['Number of rxns removed in this iter: ',num2str(numel(inactive))])    
        count=count+numel(inactive);   
        display(['Number of rxns removed so far: ',num2str(count)])
        core_inactivated=intersect(core,inactive); 
		core_removed=union(core_removed,core_inactivated);
        core=setdiff(core,inactive);
        core_active_g=setdiff(core_active_g,inactive);
        display(['Number of core rxns removed: ',num2str(numel(core_inactivated))])   
        display(['Number of rxns to remove: ',num2str(numel(rxns2remove))])  
        %continue;
    end
    % Reactions not present in any of the tissue samples based on expression
    % data form a high confidence negative set. mCADRE tries to remove them
    % even if some core reactions have to be removed, if 3 conditions are
    % met, as commented above
    partofcore=intersect(inactive,core_active_g);
    inactive=setdiff(inactive,core);
    temp_model=removeRxns(tissue_model,inactive);
    if numel(partofcore)==0&&checkModelFunction(temp_model,precursorMets)&&checkSalvage(temp_model)
        tissue_model=temp_model;
        rxns2remove=rxns2remove(~ismember(rxns2remove,inactive),1);%changed on May-09-2012
        display(['Number of rxns removed in this iter: ',num2str(numel(inactive))])    
        count=count+numel(inactive);
    else
        rxns2remove(1)=[];
        display(['Number of rxns removed in this iter: 0'])
    end
    display(['Number of rxns removed so far: ',num2str(count)])
    display(['Number of rxns to remove: ',num2str(numel(rxns2remove))])
    %For all other non-core reactions, flux through the core reactions, and
    %core model functionality must be retained.
  
end

%% Step 9:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the Tissue Model. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rxnsInTissue=tissue_model.rxns;
tissue_model_originalID=removeRxns(originalModel,setdiff(originalModel.rxns,rxnsInTissue));
%tissue_model_originalID retains the original gene id format in Recon3D

tissue_model=removeNonUsedGenes(tissue_model);
tissue_model_originalID=removeNonUsedGenes(tissue_model_originalID);
%removeNonUsedGenes removes genes no longer used in gene-reaction rules
%after the corresponding reactions are removed.
[c,ind]=ismember(tissue_model.rxns,generxns);
gRxns=tissue_model.rxns(c,:);

save TissueReconResults tissue_model tissue_model_originalID core_removed 
else
    display('genericModel failed to pass precursor metabolites test!');
end












