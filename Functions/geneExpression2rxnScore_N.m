function [geneAssociatedRxns, reactionScores] = geneExpression2rxnScore_N(expressionid, expression, parsedGPR, corrRxn)
    numGenes = numel(expressionid);
    uniqueRulesGenes = unique(parsedGPR);
    uniqueRulesGenes(1) = []; % Remove empty cell

    % Convert parsedGPR to parsedGPRmat
    parsedGPRmat = cell2mat(cellfun(@(x) str2double(x), parsedGPR, 'UniformOutput', false));

    % Create a gene-to-expression index map
    [geneIdx, expressionIdx] = ismember(parsedGPRmat, expressionid);

    % Create expressionMat using the gene-to-expression index map
    expressionMat = ones(size(parsedGPR));
    expressionMat(geneIdx) = expression(expressionIdx(geneIdx));

    % Calculate expressionCorrRxn and uniqueCorrRxn
    expressionCorrRxn = min(expressionMat, [], 2);
    uniqueCorrRxn = unique(corrRxn);

    % Compute uniqueCorrRxnScore
    reactionScores = accumarray(findgroups(corrRxn), expressionCorrRxn, [], @max);

    geneAssociatedRxns = uniqueCorrRxn;
end
