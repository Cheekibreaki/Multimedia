function QP = findCorrectQP(per_block_row_budget, bitCountPerRow)
    % Convert variables to int32 if necessary
    per_block_row_budget = int32(per_block_row_budget);
    bitCountPerRow = int32(bitCountPerRow);
    
    % Iterate through the bitCountPerRow array up to the second last element
    for idx = 1:(length(bitCountPerRow) - 1)
        val = bitCountPerRow(idx + 1);
        next_val = bitCountPerRow(idx);
        % Check if the budget falls within the current range
        if (per_block_row_budget > val && per_block_row_budget <= next_val)
            QP = idx - 1; % Adjust for 0-based indexing
            return;
        end
    end
    % Check if per_block_row_budget is greater than or equal to the last element
    if per_block_row_budget <= bitCountPerRow(end)
        QP = length(bitCountPerRow) - 1; % Adjust for 0-based indexing
     
    end
    % If no match is found, return -1
    if(per_block_row_budget > bitCountPerRow(1))
        QP = 0;
    end
end