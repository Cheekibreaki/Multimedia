function [reconstructedResiduals] = invquantization_block(quantizedResiduals, dct_blockSize, width, height, baseQP, vbs_matrix)
    [vbsRows, vbsCols] = size(vbs_matrix);
    reconstructedResiduals = zeros(size(quantizedResiduals));

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
            
            % Check if the vbs_block is all zeros or all ones
            if all(vbs_block(:) == 0)
                % Process as a large block
                % Extract the quantized residual block
                quantizedBlock = quantizedResiduals(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1);
                
                % Create the quantization matrix
                Q = createQMatrix(size(quantizedBlock), baseQP);
                
                % Inverse quantization (element-wise multiplication)
                dequantizedBlock = quantizedBlock .* Q;
                
                % Apply inverse DCT to the dequantized block
                idctBlock = idct2(dequantizedBlock);
                
                % Store the result in the reconstructed residuals
                reconstructedResiduals(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = idctBlock;
            elseif all(vbs_block(:) == 1)
                % Process each sub-block separately
                for subBlockY = 0:1
                    for subBlockX = 0:1
                        % Calculate sub-block offsets
                        subRowOffset = rowOffset + subBlockY * dct_blockSize;
                        subColOffset = colOffset + subBlockX * dct_blockSize;
                        
                        % Calculate actual sub-block size
                        actualSubBlockHeight = min(dct_blockSize, height - subRowOffset + 1);
                        actualSubBlockWidth = min(dct_blockSize, width - subColOffset + 1);
                        
                        % Extract the quantized sub-block
                        quantizedSubBlock = quantizedResiduals(subRowOffset:subRowOffset+actualSubBlockHeight-1, subColOffset:subColOffset+actualSubBlockWidth-1);
                        
                        % Create the quantization matrix
                        Q = createQMatrix(size(quantizedSubBlock), baseQP-1);
                        
                        % Inverse quantization (element-wise multiplication)
                        dequantizedSubBlock = quantizedSubBlock .* Q;
                        
                        % Apply inverse DCT to the dequantized sub-block
                        idctSubBlock = idct2(dequantizedSubBlock);
                        
                        % Store the result in the reconstructed residuals
                        reconstructedResiduals(subRowOffset:subRowOffset+actualSubBlockHeight-1, subColOffset:subColOffset+actualSubBlockWidth-1) = idctSubBlock;
                    end
                end
            else
                error('Invalid vbs_matrix block: 2x2 block must be all zeros or all ones.');
            end
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
