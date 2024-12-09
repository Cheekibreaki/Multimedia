function [quantizedResiduals] = quantization_block(residuals, dct_blockSize, width, height, baseQP, vbs_matrix)
    [vbsRows, vbsCols] = size(vbs_matrix);
    quantizedResiduals = zeros(size(residuals));

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
                % Extract the residual block
                block = residuals(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1);
                
                % Apply DCT to the block
                dctBlock = dct2(double(block));
                
                % Create the quantization matrix
                Q = createQMatrix(size(block), baseQP);
                
                % Quantization
                quantizedBlock = round(dctBlock ./ Q);
                
                % Store the quantized block
                quantizedResiduals(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = quantizedBlock;
                
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
                        
                        % Extract the residual sub-block
                        block = residuals(subRowOffset:subRowOffset+actualSubBlockHeight-1, subColOffset:subColOffset+actualSubBlockWidth-1);
                        
                        % Apply DCT to the sub-block
                        dctBlock = dct2(double(block));
                        
                        % Create the quantization matrix
                        Q = createQMatrix(size(block), baseQP);
                        
                        % Quantization
                        quantizedBlock = round(dctBlock ./ Q);
                        
                        % Store the quantized sub-block
                        quantizedResiduals(subRowOffset:subRowOffset+actualSubBlockHeight-1, subColOffset:subColOffset+actualSubBlockWidth-1) = quantizedBlock;
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
