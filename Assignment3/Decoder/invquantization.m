function [reconstructedResiduals] = invquantization(quantizedResiduals, dct_blockSize, width, height, baseQP)
    reconstructedResiduals = zeros(size(quantizedResiduals));
    for row = 1:dct_blockSize:height
        for col = 1:dct_blockSize:width
            % Extract the current block from the quantized residuals
            quantizedBlock = quantizedResiduals(row:row+dct_blockSize-1, col:col+dct_blockSize-1);
            % QP = adaptiveQP(dct_blockSize, baseQP);
            % Recreate the quantization matrix
            Q = createQMatrix(size(quantizedBlock), baseQP);
            
            % Inverse quantization (element-wise multiplication)
            dequantizedBlock = quantizedBlock .* Q;
            

            % if rounded, will affect quality of I frame !!!!

            

            % Apply inverse DCT to the dequantized block
            idctBlock = idct2(dequantizedBlock);
            
            % Store the result in the reconstructed residuals
            reconstructedResiduals(row:row+dct_blockSize-1, col:col+dct_blockSize-1) = idctBlock;
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