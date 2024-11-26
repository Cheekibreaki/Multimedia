function test_quantization_consistency()
    % Parameters
    dct_blockSize = 8;
    width = 16;
    height = 16;
    baseQP = 5;

    % Generate synthetic residuals
    residuals = randi([-255, 255], height, width);

    % Case 1: Quantization without vbs_matrix
    quantizedResiduals_no_vbs = quantization(residuals, dct_blockSize, width, height, baseQP-1);

    % Case 2: Quantization with vbs_matrix filled with ones
    vbs_matrix_rows = height / dct_blockSize;
    vbs_matrix_cols = width / dct_blockSize;
    vbs_matrix = ones(vbs_matrix_rows, vbs_matrix_cols);

    % Modify the quantization function to use baseQP for sub-blocks
    % (Assuming you have updated the function accordingly)

    quantizedResiduals_with_vbs = quantization(residuals, dct_blockSize, width, height, baseQP, vbs_matrix);

    % Compare the outputs
    difference = quantizedResiduals_no_vbs - quantizedResiduals_with_vbs;
    max_difference = max(abs(difference(:)));

    % Display the results
    if max_difference == 0
        disp('The quantized residuals are identical for both cases.');
    else
        disp('The quantized residuals differ between the two cases.');
        fprintf('Maximum absolute difference: %d\n', max_difference);
    end
end


test_quantization_consistency()






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


function [quantizedResiduals] = quantization(residuals, dct_blockSize, width, height, baseQP, vbs_matrix)
    quantizedResiduals = zeros(size(residuals));
    if exist('vbs_matrix', 'var') && ~isempty(vbs_matrix)
    [vbsRows, vbsCols] = size(vbs_matrix);
    
    
    % Iterate over the vbs_matrix in steps of 2 to process 2x2 blocks
    for blockY = 1:2:vbsRows
        for blockX = 1:2:vbsCols
            % Extract the 2x2 block from vbs_matrix
            vbs_block = vbs_matrix(blockY:blockY+1, blockX:blockX+1);

            % Determine the starting position in the residuals
            rowOffset = (blockY - 1) * dct_blockSize + 1;
            colOffset = (blockX - 1) * dct_blockSize + 1;

            % Calculate the actual block size, accounting for image boundaries
            actualBlockHeight = min(2 * dct_blockSize, height - rowOffset + 1);
            actualBlockWidth = min(2 * dct_blockSize, width - colOffset + 1);

            % Check if the vbs_block is all ones or contains zeros
            if all(vbs_block(:) == 0)
                % Process as a large block
                % Extract the residual block
                block = residuals(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1);

                % Apply DCT to the block
                dctBlock = dct2(double(block));

                % Create the quantization matrix
                Q = createQMatrix(size(dctBlock), baseQP);

                % Quantize the DCT coefficients
                quantizedBlock = round(dctBlock ./ Q);

                % Store the quantized block
                quantizedResiduals(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = quantizedBlock;
            else
                % Process each sub-block separately
                for subBlockY = 0:1
                    for subBlockX = 0:1
                        % Calculate sub-block offsets
                        subRowOffset = rowOffset + subBlockY * dct_blockSize;
                        subColOffset = colOffset + subBlockX * dct_blockSize;

                        % Check if within image boundaries
                        if subRowOffset > height || subColOffset > width
                            continue;
                        end

                        % Calculate actual sub-block size
                        actualSubBlockHeight = min(dct_blockSize, height - subRowOffset + 1);
                        actualSubBlockWidth = min(dct_blockSize, width - subColOffset + 1);

                        % Extract the residual sub-block
                        block = residuals(subRowOffset:subRowOffset+actualSubBlockHeight-1, subColOffset:subColOffset+actualSubBlockWidth-1);

                        % Apply DCT to the sub-block
                        dctBlock = dct2(double(block));

                        % Create the quantization matrix
                        Q = createQMatrix(size(dctBlock), baseQP -1);

                        % Quantize the DCT coefficients
                        quantizedBlock = round(dctBlock ./ Q);

                        % Store the quantized sub-block
                        quantizedResiduals(subRowOffset:subRowOffset+actualSubBlockHeight-1, subColOffset:subColOffset+actualSubBlockWidth-1) = quantizedBlock;
                    end
                end
            end
        end
    end
    else
        quantizedResiduals = zeros(size(residuals));
        for row = 1:dct_blockSize:height
            for col = 1:dct_blockSize:width
                block = residuals(row:row+dct_blockSize-1, col:col+dct_blockSize-1);
                dctBlock = dct2(double(block));

                Q = createQMatrix(size(block), baseQP);
                quantizedBlock = round(dctBlock ./ Q);
                quantizedResiduals(row:row+dct_blockSize-1, col:col+dct_blockSize-1) = quantizedBlock;
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