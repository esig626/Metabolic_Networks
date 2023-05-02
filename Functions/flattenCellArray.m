function flattened = flattenCellArray(input)
    if ~iscell(input)
        flattened = input;
    else
        flattened = {};
        for i = 1:numel(input)
            temp = flattenCellArray(input{i});
            flattened = [flattened, temp];
        end
    end
end