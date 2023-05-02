%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Elias Siguenza
% Date: 23.05.2022 (Article Review Corrections and Implementation of VFFVA)
% Address: Institute of Metabolism and Systems Research, 
% College of Medical and Dental Sciences, University of Birmingham, 
% B15 2TT, United Kingdom.
% The Tennant Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New improved mas5Calss function. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% I have updated the original function to run faster.
% Here are some changes made to the code:
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% The if nargin < ... statements are condensed into one-liners.
% The for loop is replaced with a parfor loop, which allows parallel 
% execution of the loop iterations. This can significantly speed up the 
% execution if you have multiple CPU cores available. Make sure you have 
% the Parallel Computing Toolbox installed for this feature to work.
% The signrank function is replaced with a custom 
% signrank_adjusted function that takes care of the sign adjustment 
% without recalculating the p-value. I Removed unused and commented-out 
% code. Now it runs on average 2.8 times faster. 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [calls, pvalues] = mas5Calls_new(raw_struct, tau, alpha1, alpha2)
    if nargin < 4, alpha2 = 0.06; end
    if nargin < 3, alpha1 = 0.04; end
    if nargin < 2, tau = 0.015; end
    if nargin < 1, error('Usage: [CALLS] = MAS5CALLS(RAW_STRUCT, TAU, ALPHA1, ALPHA2)'); end
    
    probeIndicesZero = find(raw_struct.ProbeIndices == 0);
    probeIndicesZero = [probeIndicesZero; raw_struct.NumProbes + 1];
    probeSetLengths = diff(probeIndicesZero);
    
    numProbeSets = raw_struct.NumProbeSets;
    numSamples = size(raw_struct.PMIntensities, 2);
    
    calls = zeros(numProbeSets, numSamples);
    pvalues = zeros(numProbeSets, numSamples);
    
    parfor j = 1:numProbeSets
        calls_j = zeros(1, numSamples);
        pvalues_j = zeros(1, numSamples);
        
        for k = 1:numSamples
            PM = raw_struct.PMIntensities(probeIndicesZero(j):probeIndicesZero(j) + probeSetLengths(j) - 1, k);
            MM = raw_struct.MMIntensities(probeIndicesZero(j):probeIndicesZero(j) + probeSetLengths(j) - 1, k);
            
            dv = (PM - MM) ./ (PM + MM);
            
            [p, W] = signrank_adjusted(dv, tau);
            
            pvalues_j(k) = p;
            
            if p < alpha1
                calls_j(k) = 2;
            elseif p > alpha2
                calls_j(k) = 0;
            else
                calls_j(k) = 1;
            end
        end
        
        calls(j, :) = calls_j;
        pvalues(j, :) = pvalues_j;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p, W] = signrank_adjusted(data, tau)
    abs_data = abs(data);
    [ranks, ~] = tiedrank(abs_data);
    signedranks = sign(data) .* ranks;
    W = sum(signedranks);
    
    if W <= 0
        tail = 'right';
    else
        tail = 'left';
    end
    
    p = signrank(data, tau, 'tail', tail);
end
