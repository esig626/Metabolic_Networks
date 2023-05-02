function [parsedGPR,corrRxn,corrRxn_iter] = extractGPRs_N(model)

warning off all

parsedGPR = [];
corrRxn = [];
cnt = 1;

for i = 1:length(model.rxns)
    if length(model.grRules{i}) > 1
        % Parsing each reactions gpr
        [parsing{1,1},parsing{2,1}] = strtok(model.grRules{i},'or');
        for j = 2:1000
            [parsing{j,1},parsing{j+1,1}] = strtok(parsing{j,1},'or');
            if isempty(parsing{j+1,1})==1
                break
            end
        end
        
        for j = 1:length(parsing)
            for k = 1:1000
                [parsing{j,k},parsing{j,k+1}] = strtok(parsing{j,k},'and');
                if isempty(parsing{j,k+1})==1
                    break
                end
            end
        end
        
        for j = 1:size(parsing,1)
            for k = 1:size(parsing,2)
                parsing{j,k} = strrep(parsing{j,k},'(','');
                parsing{j,k} = strrep(parsing{j,k},')','');
                parsing{j,k} = strrep(parsing{j,k},' ','');
            end
        end
        
        for j = 1:size(parsing,1)-1
            newparsing(j,:) = parsing(j,1:length(parsing(j,:))-1);
        end
        
        parsing = newparsing;
        
     
        for j = 1:size(parsing,1)
            for k = 1:size(parsing,2)
                if length(parsing{j,k}) == 0
                    parsing{j,k} = '';                    
                end
            end
        end
        
              
        num = size(parsing,1);
        for j = 1:num
            sizeP = length(parsing(j,:));
            if sizeP > size(parsedGPR,2)
                for k = 1:size(parsedGPR,1)
                    parsedGPR{k,sizeP} = {''};
                end
            end
            
            for l = 1:sizeP          
            parsedGPR{cnt,l} = parsing(j,l);
            end           
            cnt = cnt+1;
        end
        
        for j = 1:num
            corrRxn = [corrRxn;model.rxns(i)];
        end
        
        clear parsing newparsing
        
    end
end

% % Find empty rows
empty_rows = cellfun(@(x) strcmp(x, ''), parsedGPR(:, 1));
corrRxn(empty_rows, :) = [];
parsedGPR(empty_rows, :) = [];
% Remove empty rows
for i = 1:size(parsedGPR,1)
    for j = 1:size(parsedGPR,2)
        if isempty(parsedGPR{i,j}) == 1
            parsedGPR{i,j} = {''};
        end
    end
end



for i = 1:size(parsedGPR,1)
    for j= 1:size(parsedGPR,2)
        if isempty(parsedGPR{i,j})
            parsedGPR2(i,j) = {''};
        else
            parsedGPR2(i,j) = parsedGPR(i,j);
        end
    end
end

parsedGPR = parsedGPR2;

end
