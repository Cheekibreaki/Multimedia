function [quantizedResiduals] = quantization(residuals, dct_blockSize,width,height,baseQP)
    % if(RCflag == false) 
        quantizedResiduals = zeros(size(residuals));
        for row = 1:dct_blockSize:height
            for col = 1:dct_blockSize:width
                block = residuals(row:row+dct_blockSize-1, col:col+dct_blockSize-1);
                dctBlock = dct2(double(block));
        
                % Quantization
                % QP = adaptiveQP(dct_blockSize, baseQP);
                Q = createQMatrix(size(block), baseQP);
                quantizedBlock = round(dctBlock ./ Q);
                quantizedResiduals(row:row+dct_blockSize-1, col:col+dct_blockSize-1) = quantizedBlock;
            end
        end
    % else
    %     curr_per_block_row_budget = per_block_row_budget;
    %     quantizedResiduals = zeros(size(residuals));
    %     for row = 1:dct_blockSize:height
    %         QP = findCorrectQP(curr_per_block_row_budget, bitCountPerRow);
    %         for col = 1:dct_blockSize:width
    %             block = residuals(row:row+dct_blockSize-1, col:col+dct_blockSize-1);
    %             dctBlock = dct2(double(block));
    % 
    %             % Quantization
    %             % QP = adaptiveQP(dct_blockSize, baseQP);
    %             Q = createQMatrix(size(block), baseQP);
    %             quantizedBlock = round(dctBlock ./ Q);
    %             quantizedResiduals(row:row+dct_blockSize-1, col:col+dct_blockSize-1) = quantizedBlock;
    %         end
    %         find the size 
    %     end
    % end

end
function QP = findCorrectQP(per_block_row_budget, bitCountPerRow)
    % Iterate through the bitCountPerRow array
    for idx = 1:(length(bitCountPerRow))
        % Check if the budget falls within the current range
        if idx == length(bitCountPerRow)
            QP = idx - 1;
            return;
        end
        if per_block_row_budget >= bitCountPerRow(idx) && per_block_row_budget < bitCountPerRow(idx + 1)
            QP = idx - 1; % Return the corresponding QP value (adjust for 0-based indexing)
            return;
        end
    end
    % If no match is found, return -1 or an error message
    QP = -1; % Indicates that no valid QP was found
end

function Q = createQMatrix(blockSize, QP)
    if numel(blockSize) > 1
        rows = blockSize(1);
        cols = blockSize(2);
    else
        rows = blockSize;
        cols = blockSize;
    end
    
    Q = zeros(rows, cols);
    for x = 1:rows
        for y = 1:cols
            if (x + y <= rows)
                Q(x,y) = 2^QP;
            elseif (x + y == rows + 1)
                Q(x,y) = 2^(QP+1);
            else
                Q(x,y) = 2^(QP+2);
            end
        end
    end
end

% function QP = adaptiveQP(i, baseQP)
%     % Check if QP is in the valid range
% 
% 
%     % Initialize the quantization matrix Q
%     QP = zeros(i, i);
%     for row = 1:i
%         for col = 1:i
%             QP(row, col) = 2 ^ (baseQP + floor((row + col - 2) / i));
%         end
%     end
% end