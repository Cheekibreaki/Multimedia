function [quantizedResiduals] = quantization(residuals, dct_blockSize,width,height,baseQP)
    quantizedResiduals = zeros(size(residuals));
    for row = 1:dct_blockSize:height
        for col = 1:dct_blockSize:width
            block = residuals(row:row+dct_blockSize-1, col:col+dct_blockSize-1);
            dctBlock = dct2(double(block));
    
            % Quantization
            % QP = adaptiveQP(dct_blockSize, baseQP);
            Q = createQMatrix(size(block), baseQP);
            quantizedBlock = dctBlock ./ Q;
            quantizedResiduals(row:row+dct_blockSize-1, col:col+dct_blockSize-1) = quantizedBlock;
        end
    end
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