function [quantizedResiduals] = quantization(residuals, dct_blockSize,width,height,QP)
    quantizedResiduals = zeros(size(residuals));
    for row = 1:dct_blockSize:height
        for col = 1:dct_blockSize:width
            block = residuals(row:row+dct_blockSize-1, col:col+dct_blockSize-1);
            dctBlock = dct2(double(block));
    
            % Quantization
            Q = createQMatrix(size(block), QP);
            quantizedBlock = (dctBlock ./ Q);
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

% function adaptiveQP = getAdaptiveQP(block, baseQP)
%     blockVariance = var(double(block(:)));
%     if blockVariance > 1000
%         adaptiveQP = max(0, baseQP - 2);  % Reduce QP for high-variance blocks
%     elseif blockVariance < 100
%         adaptiveQP = min(51, baseQP + 2);  % Increase QP for low-variance blocks
%     else
%         adaptiveQP = baseQP;
%     end
% end