function [inactiveRxns,time,numOpts] = checkModelConsistency_N(model,r,C,fast,path)

% INPUT
%  model                COBRA model structure
%  modelmps             Mathematical Programing Script version of the model
%  r                    name of reaction to be removed (for model pruning)
%  C                    list of core reaction names (for model pruning)
%  path                 Path to save the .mps file for VFFVA 




% When using the function independently, no additional inputs need to be
% specified besides the model.
% nargin = number of arguments in.

if nargin < 2
    r = [];
    C = {};
end
if nargin < 4
    fast = 0;
end
if numel(r)
   % Remove reaction r from the model
    disp(['Model Length:',num2str(length(model.rxns))])
    model = removeRxns(model,r); 
    disp(['New Length:', num2str(length(model.rxns))])
    
end

if fast==0
    model.c(logical(model.c)) = 0;
    
else
    % Change the processing and paralellizing parameters of the VFFVA
    % routine here. 
    nCores=3;
    nThreads=2;
    ex='';
end

inactiveRxns = r;

numOpts = 0;
t0 = clock;

% First check whether any core reactions are blocked by the removal of r.
rxnList = C;

% Checking for metabolite dead ends is accomplished entirely by matrix
% operations, and is therefore very fast in Matlab. If any core reaction
% contains a dead-end metabolite, the reaction itself will be a dead end.
% This check potentially avoids sequential optimizations, as the function
% can exit if any blocked core reactions are detected.
deadEndMets = detectDeadEnds_fast(model);
deadEnd = sum(full(model.S(deadEndMets,:)~=0),1)>0;
deadEndRxns = model.rxns(deadEnd);
% Added legacy option because intersect has changed in latest MATLAB
% version.
deadEnd_C = intersect(rxnList,deadEndRxns,'legacy');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% This works when removing reactions at the last step of
%%%%%%%%%%%%%%%%%%%% the Main File.
% The original mCadre sets objectives of the
% optimisation to zero. This effectively means that you are removing any 
% objective function from the model, and the subsequent FVA calculation 
% will find the maximum and minimum flux values for each reaction while 
% maintaining a feasible solution that satisfies the constraints but 
% without optimizing towards any specific objective.

% % VFFVA does not accept a zero objective function (but does so inside the C
% % file) to get around that, and without changing the C code, we can set a
% % redundant reaction (TESTED!), which does not affect the optimisation 
% % of the model: 
model = addReaction(model,'FVA_zero','reactionFormula','A -> B + 2 C');
model.c(:)=0; model.c(end)=1;
% We now recompute the model as an .mps (Mathematical Programming System)
cd(path)
modelmps = model;
output_filename = 'model';
convertProblem(modelmps, 0, output_filename);
modelmps = 'model.mps';
cd ..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(deadEnd_C)
    inactiveRxns = union(deadEnd_C,inactiveRxns);
elseif fast
    disp('Hold on to your seat, I am initialising the system ...')
    %Set install folder of MPI
    setenv('PATH', ['/opt/local/bin:', getenv('PATH')]);
    [status, cmdout] = system('mpiexec --version');
    %Set schedule and chunk size parameters
    setenv('OMP_SCHEDUELE', ['dynamic' ',' num2str(50)])
    % Get the .mps formatted model 
    cd(path); 
    clc
    % Run VFFVA!
    disp('Running Very Fast Flux Variability Analysis!')
    [status,cmdout]=system(['mpirun -np ' num2str(nCores)...
    ' --bind-to ' num2str('none') ' -x OMP_NUM_THREADS=' num2str(nThreads)...
    ' ./veryfastFVA ' modelmps ' ' num2str(90)  ' ' num2str(0) ' ' ex]);
    clc
    disp('I am done! I will now fetch your results, hang tight...')
    resultFile=[modelmps(1:end-4) 'output.csv'];
    results=readtable(resultFile);
    optMin=results.minFlux; optMax=results.maxFlux;
    fastInactive = (abs(optMax) < 1e-6) & (abs(optMin) < 1e-6);
    fastInactiveRxns = model.rxns(fastInactive);
    inactiveRxns = union(inactiveRxns,fastInactiveRxns);
    clc; disp('All done!')
    cd ..
    


% Othwerise, if you dont have fastFVA or Very Fast FVA the next 
% step follows the heuristic speed-up to FVA proposed by
% Jerby et al. to scan through core reactions only; if any core reactions
% are blocked, the function can exit.
else
    disp('Checking core reactions...')
    changeCobraSolver('glpk');
    clc
    while numel(rxnList)
        numRxnList = numel(rxnList);
        model = changeObjective(model,rxnList);

        % Maximize all
        FBAsolution = optimizeCbModel(model,'max');
        numOpts = numOpts + 1;
        optMax = FBAsolution.x;
        active = (abs(optMax) >= 1e-6);
        activeRxns = model.rxns(active);
        rxnList = setdiff(rxnList,activeRxns);

        % Minimize all
        if ~numel(activeRxns)
            FBAsolution = optimizeCbModel(model,'min');
            numOpts = numOpts + 1;
            optMin = FBAsolution.x;
            active = (abs(optMin) >= 1e-6);
            activeRxns = model.rxns(active);
            rxnList = setdiff(rxnList,activeRxns);
        end

        numRemoved = numRxnList - numel(rxnList);

        if ~numRemoved
            randInd = randperm(numel(rxnList));
            i = rxnList(randInd(1));
            model = changeObjective(model,i);

            % Maximize reaction i
            FBAsolution = optimizeCbModel(model,'max');
            numOpts = numOpts + 1;
            optMax = FBAsolution.f;

            % Minimize reaction i (if reversible)
            lb_i = model.lb(strmatch(i,model.rxns,'exact'));
            if lb_i < 0 && optMax < 1e-6
                FBAsolution = optimizeCbModel(model,'min');
                numOpts = numOpts + 1;
                optMin = FBAsolution.f;
            else optMin = 0;
            end

            if (abs(optMax) < 1e-6) && (abs(optMin) < 1e-6)
                inactiveRxns = union(inactiveRxns,i);
                break;
            end

            rxnList = setdiff(rxnList,i);
        end
    end
end

% If all core actions remain active, r can be removed from the model. This
% step checks all other non-core reactions to determine whether anything
% else can be concurrently removed. Note: core reactions are included in
% the objective function, as this seems to speed up the algorithm, but are
% not individuall maximized or minimized (since this is redundant witht he
% above step).
if numel(inactiveRxns) <= 1 && fast == 0  
    disp('Checking non-core reactions...')
    rxnList = model.rxns;
    inactiveRxns = union(inactiveRxns,deadEndRxns);
    rxnList = setdiff(rxnList,inactiveRxns);

    while numel(rxnList)        
        numRxnList = numel(rxnList);
        model = changeObjective(model,rxnList);

        % Maximize all
        FBAsolution = optimizeCbModel(model,'max');
        numOpts = numOpts + 1;
        optMax = FBAsolution.x;
        active = (abs(optMax) >= 1e-6);
        activeRxns = model.rxns(active);
        rxnList = setdiff(rxnList,activeRxns);

        % Minimize all
        if ~numel(activeRxns)
            FBAsolution = optimizeCbModel(model,'min');
            numOpts = numOpts + 1;
            optMin = FBAsolution.x;
            active = (abs(optMin) >= 1e-6);
            activeRxns = model.rxns(active);
            rxnList = setdiff(rxnList,activeRxns);
        end

        numRemoved = numRxnList - numel(rxnList);

        if ~numRemoved
            randInd = randperm(numel(rxnList));
            i = rxnList(randInd(1));
            
            if isempty(C)
                C={};
            end
            
            if ~ismember(i,C)
                model = changeObjective(model,i);

                % Maximize reaction i
                FBAsolution = optimizeCbModel(model,'max');
                numOpts = numOpts + 1;
                optMax = FBAsolution.f;

                % Minimize reaction i (if reversible)
                lb_i = model.lb(strmatch(i,model.rxns,'exact'));
                if lb_i < 0 && optMax < 1e-6
                    FBAsolution = optimizeCbModel(model,'min');
                    numOpts = numOpts + 1;
                    optMin = FBAsolution.f;
                else optMin = 0;
                end

                if (abs(optMax) < 1e-6) && (abs(optMin) < 1e-6)%4.462525
                    inactiveRxns = union(inactiveRxns,i);
                end
            end

            rxnList = setdiff(rxnList,i);
            if ~numel(setdiff(rxnList,C))
                break;
            end
        end
    end
end

time = etime(clock,t0);
display(['checkModelConsistency time: ',num2str(time,'%1.2f'),' s'])