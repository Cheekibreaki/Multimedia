

 
function [motionVectors, avgMAE,vbs_matrix] = vbs_motionEstimation(currentFrame, referenceFrames, blockSize, searchRange, dct_blockSize,QP,lambda)
    % Motion Estimation function that processes blocks in raster order.
    % Calls findBestMatch to find the best matching block in the reference frame.
    %
    % Parameters:
    %   currentFrame - The Y component of the current frame (grayscale image)
    %   referenceFrames - The Y component of the reference frames (grayscale images)
    %   blockSize - The size of the blocks (e.g., 8x8)
    %   searchRange - The maximum number of pixels to search in any direction for motion
    %   VBS_matrix - Matrix containing the VBS information (2 for single large block, 1 for four sub-blocks)
    %
    % Returns:
    %   motionVectors - An array of motion vectors for each block
    %   avgMAE - The average Mean Absolute Error for the entire frame

    % Get the dimensions of the frame
    [height, width] = size(currentFrame); % Example: 288x352

    % Initialize the motion vector array (for storing motion vectors for each block)
    numBlocksX = width / blockSize;
    numBlocksY = height / blockSize;
    motionVectors = zeros(numBlocksY, numBlocksX, 3);  % Stores motion vector for each block.
    vbs_matrix = -1 * ones(numBlocksY, numBlocksX);
    totalMAE = 0;  % To keep track of the total MAE across all blocks


    for blockY = 1:2:numBlocksY
        for blockX = 1:2:numBlocksX

            % Single large block: size is 2 * blockSize (i.e., twice the size in both dimensions)
            currentBlockSize = blockSize * 2;
            rowOffset = (blockY - 1) * blockSize + 1;
            colOffset = (blockX - 1) * blockSize + 1;

            % Extract the current block from the current frame
            currentBlock = currentFrame(rowOffset:rowOffset + currentBlockSize - 1, ...
                                        colOffset:colOffset + currentBlockSize - 1);

            % Compute motion estimation for large block and split blocks
            [motionVector_block_large, total_minMAE_large] = compute_motionVector_block(currentBlock, currentBlockSize, referenceFrames, rowOffset, colOffset, searchRange, true);
            [motionVector_block_split, total_minMAE_split] = compute_motionVector_block(currentBlock, currentBlockSize, referenceFrames, rowOffset, colOffset, searchRange, false);

           
            predictedFrame_block_large = compute_predictedFrame_block(referenceFrames, motionVector_block_large, rowOffset, colOffset, blockSize);
            predictedFrame_block_split = compute_predictedFrame_block(referenceFrames, motionVector_block_split, rowOffset, colOffset, blockSize);
            
            Residuals_block_large = double(currentBlock) - double(predictedFrame_block_large);
            Residuals_block_split = double(currentBlock) - double(predictedFrame_block_split);

            quantizedResiduals_large = quantization(Residuals_block_large, dct_blockSize,currentBlockSize,currentBlockSize,QP);
            
             % Split the 16x16 residue into four 8x8 blocks
            blockSizeSmall = blockSize; % blockSize here is 8
            Residuals_split_1 = Residuals_block_split(1:blockSizeSmall, 1:blockSizeSmall);
            Residuals_split_2 = Residuals_block_split(1:blockSizeSmall, blockSizeSmall+1:end);
            Residuals_split_3 = Residuals_block_split(blockSizeSmall+1:end, 1:blockSizeSmall);
            Residuals_split_4 = Residuals_block_split(blockSizeSmall+1:end, blockSizeSmall+1:end);

            % Quantize each of the four 8x8 blocks separately
            quantizedResiduals_1 = quantization(Residuals_split_1, dct_blockSize, blockSizeSmall, blockSizeSmall, QP-1);
            quantizedResiduals_2 = quantization(Residuals_split_2, dct_blockSize, blockSizeSmall, blockSizeSmall, QP-1);
            quantizedResiduals_3 = quantization(Residuals_split_3, dct_blockSize, blockSizeSmall, blockSizeSmall, QP-1);
            quantizedResiduals_4 = quantization(Residuals_split_4, dct_blockSize, blockSizeSmall, blockSizeSmall, QP-1);

            % Combine the quantized 8x8 blocks back into a 16x16 block
            quantizedResiduals_split = [quantizedResiduals_1, quantizedResiduals_2;
                                       quantizedResiduals_3, quantizedResiduals_4];


            compresiduals_large = invquantization(quantizedResiduals_large, dct_blockSize,currentBlockSize,currentBlockSize,QP);

            
            % Compute the inverse quantized residuals for the split blocks
            compresiduals_1 = invquantization(quantizedResiduals_1, dct_blockSize, blockSizeSmall, blockSizeSmall, QP-1);
            compresiduals_2 = invquantization(quantizedResiduals_2, dct_blockSize, blockSizeSmall, blockSizeSmall, QP-1);
            compresiduals_3 = invquantization(quantizedResiduals_3, dct_blockSize, blockSizeSmall, blockSizeSmall, QP-1);
            compresiduals_4 = invquantization(quantizedResiduals_4, dct_blockSize, blockSizeSmall, blockSizeSmall, QP-1);

            % Combine the inverse quantized 8x8 blocks back into a 16x16 block
            compresiduals_split = [compresiduals_1, compresiduals_2;
                                   compresiduals_3, compresiduals_4];

            
            reconstructed_split = predictedFrame_block_split + compresiduals_split;
            reconstructed_large = predictedFrame_block_large + compresiduals_large;
            

            % Compute SAD for reconstructed_split
            SAD_split = sum(sum(abs(double(currentBlock) - double(reconstructed_split))));
            
            % Compute SAD for reconstructed_large
            SAD_large = sum(sum(abs(double(currentBlock) - double(reconstructed_large))));
            
            MDiffMV_large = diffEncoding(motionVector_block_large,'mv');
            MDiffMV_split = diffEncoding(motionVector_block_split,'mv');

            [encodedMDiff_large,nonimporatant1,encodedResidues_large] = entropyEncode(false, MDiffMV_large, [], quantizedResiduals_large);
            [encodedMDiff_split,nonimporatant1,encodedResidues_split] = entropyEncode(false, MDiffMV_split, [], quantizedResiduals_split);
            
           
            % Calculate rate (R) for large and split blocks
            R_large = (numel(encodedMDiff_large) + numel(encodedResidues_large)) / ...
                      (numel(encodedMDiff_split) + numel(encodedResidues_split) + numel(encodedMDiff_large) + numel(encodedResidues_large));

            R_split = (numel(encodedMDiff_split) + numel(encodedResidues_split)) / ...
                      (numel(encodedMDiff_split) + numel(encodedResidues_split) + numel(encodedMDiff_large) + numel(encodedResidues_large));

            D_large = SAD_large/(SAD_large+SAD_split)
            D_split = SAD_split/(SAD_large+SAD_split)
            
            J_large = D_large + lambda*R_large
            J_split = D_split + lambda*R_split
            
            % Choose the motion vector block with the lower MAE
            if J_large < J_split
                motionVector_block = motionVector_block_large;
                total_minMAE = total_minMAE_large;
                vbs_matrix(blockY:blockY+1, blockX:blockX+1) = 1;
            else
                motionVector_block = motionVector_block_split;
                total_minMAE = total_minMAE_split;
                vbs_matrix(blockY:blockY+1, blockX:blockX+1) = 0;
            end

            % Store the motion vectors in the corresponding positions
            motionVectors(blockY:blockY+1, blockX:blockX+1, :) = motionVector_block;

            % Update the total MAE
            totalMAE = totalMAE + total_minMAE;
        end
    end

    % Calculate the average MAE across all blocks
    avgMAE = totalMAE / (numBlocksX * numBlocksY /4);
end





function predictedBlock = compute_predictedFrame_block(referenceFrames, motionVector_block, rowOffset, colOffset, blockSize)
    % This function computes the predicted block from the reference frames based on motion vectors.
    %
    % Parameters:
    %   referenceFrames - Cell array of reference frames (grayscale images)
    %   motionVector_block - Motion vector for each block or sub-block (2x2x3 matrix)
    %   rowOffset - The starting row position of the large block in the frame
    %   colOffset - The starting column position of the large block in the frame
    %   blockSize - Size of the sub-block (half of the large block)
    %
    % Returns:
    %   predictedBlock - The predicted block based on the motion vectors and reference frames

    % Initialize predicted block to zeros (2 * blockSize x 2 * blockSize)
    predictedBlock = zeros(2 * blockSize, 2 * blockSize);

    % Iterate over each sub-block (2x2)
    for blockY = 1:2
        for blockX = 1:2
            % Extract motion vector components for the current sub-block
            mvY = motionVector_block(blockY, blockX, 1);
            mvX = motionVector_block(blockY, blockX, 2);
            refIdx = motionVector_block(blockY, blockX, 3) + 1; % Reference frame index

            % Calculate the position of the current sub-block in the reference frame
            subBlockRowOffset = rowOffset + (blockY - 1) * blockSize;
            subBlockColOffset = colOffset + (blockX - 1) * blockSize;
            
            % Get the corresponding reference frame
            referenceFrame = referenceFrames{refIdx};
            
            % Calculate the position of the reference sub-block based on the motion vector
            refRowOffset = subBlockRowOffset + mvY;
            refColOffset = subBlockColOffset + mvX;
            
            % Ensure the reference sub-block is within the bounds of the reference frame
            [height, width] = size(referenceFrame);
            refRowOffset = max(1, min(refRowOffset, height - blockSize + 1));
            refColOffset = max(1, min(refColOffset, width - blockSize + 1));

            % Extract the reference sub-block
            referenceBlock = referenceFrame(refRowOffset:refRowOffset + blockSize - 1, ...
                                           refColOffset:refColOffset + blockSize - 1);

            % Place the reference sub-block into the predicted block
            predictedBlock((blockY - 1) * blockSize + 1:blockY * blockSize, ...
                           (blockX - 1) * blockSize + 1:blockX * blockSize) = referenceBlock;
        end
    end
end




function [motionVector_block,total_minMAE] = compute_motionVector_block(currentBlock, currentBlockSize, referenceFrames, rowOffset, colOffset, searchRange, isLarge)
    
    % Initialize variables for the best match
    bestVector = [0, 0];
    bestRefIdx = 0;
    minMAE = inf;
    bestL1Norm = inf;
    total_minMAE = 0;
    % Split current block into 4 pieces and find the best match for each
    motionVector_block = zeros(2, 2, 3);
    subBlockSize = currentBlockSize / 2;

    if isLarge

        % Check all reference frames to find the best match
        for refIdx = 1:length(referenceFrames)
            referenceFrame = referenceFrames{refIdx};
    
            % Find the best match within the current reference frame
            [vector, mae, L1Norm] = findBestMatch(currentBlock, referenceFrame, rowOffset, colOffset, currentBlockSize, searchRange);
    
            % Update if a better match is found
            if mae < minMAE || (mae == minMAE && L1Norm < bestL1Norm)
                minMAE = mae;
                bestVector = vector;
                bestL1Norm = L1Norm;
                bestRefIdx = refIdx;
            end
        end
    
        % Store the motion vector for this large block in all corresponding sub-blocks
        motionVector_block = zeros(2, 2, 3);
        for dy = 0:1
            for dx = 0:1
                currentBlockY = dy + 1;
                currentBlockX = dx + 1;
    
                motionVector_block(currentBlockY, currentBlockX, 1) = bestVector(1);
                motionVector_block(currentBlockY, currentBlockX, 2) = bestVector(2);
                motionVector_block(currentBlockY, currentBlockX, 3) = bestRefIdx - 1;
            end
        end
        total_minMAE = minMAE;
    else


        for dy = 0:1
            for dx = 0:1
                % Extract the sub-block
                subBlockY = dy * subBlockSize + 1;
                subBlockX = dx * subBlockSize + 1;
                subBlock = currentBlock(subBlockY:subBlockY + subBlockSize - 1, subBlockX:subBlockX + subBlockSize - 1);
                
                % Initialize variables for the best match of the sub-block
                bestVector = [0, 0];
                bestRefIdx = 0;
                minMAE = inf;
                bestL1Norm = inf;

                % Check all reference frames to find the best match for the sub-block
                for refIdx = 1:length(referenceFrames)
                    referenceFrame = referenceFrames{refIdx};
                    
                    % Find the best match for the sub-block within the current reference frame
                    [vector, mae, L1Norm] = findBestMatch(subBlock, referenceFrame, rowOffset + subBlockY - 1, colOffset + subBlockX - 1, subBlockSize, searchRange);
                    
                    % Update if a better match is found
                    if mae < minMAE || (mae == minMAE && L1Norm < bestL1Norm)
                        minMAE = mae;
                        bestVector = vector;
                        bestL1Norm = L1Norm;
                        bestRefIdx = refIdx;
                    end
                end
                
                % Store the motion vector for the sub-block
                currentBlockY = dy + 1;
                currentBlockX = dx + 1;
                motionVector_block(currentBlockY, currentBlockX, 1) = bestVector(1);
                motionVector_block(currentBlockY, currentBlockX, 2) = bestVector(2);
                motionVector_block(currentBlockY, currentBlockX, 3) = bestRefIdx - 1;

                total_minMAE = total_minMAE + minMAE;
            end
        end
        total_minMAE = total_minMAE/4;
    end
end




function [bestVector, minMAE, bestL1Norm] = findBestMatch(currentBlock, referenceFrame, row, col, blockSize, searchRange)
    % Find the best matching block using Mean Absolute Error (MAE).
    %
    % Parameters:
    %   currentBlock - The block from the current frame to match
    %   referenceFrame - The reference frame 
    %   row, col - The starting position of the block in the current frame
    %   blockSize - The size of the block
    %   searchRange - The maximum pixel offset to search in the reference frame
    %
    % Returns:
    %   bestVector - The best motion vector (dy, dx) for the current block
    %   minMAE - The minimum MAE for the best matching block


    % Initialization
    minMAE = inf; 
    bestVector = [0, 0]; 
    bestL1Norm =inf;

    % Reference Frame boundry
    [heightBoundry,widthBoundry] = size(referenceFrame);

    % Iterate over the search range to find the best matching block in the reference frame
    for yOffset = -searchRange:searchRange
        for xOffset = -searchRange:searchRange
            % Determine the reference block's starting coordinates
            refRow = row + yOffset;
            refCol = col + xOffset;

            % Check if the reference block is within bounds
            if refRow > 0 && refRow + blockSize - 1 <= heightBoundry && refCol > 0 && refCol + blockSize - 1 <= widthBoundry
                % Extract the reference block
                refBlock = referenceFrame(refRow:refRow + blockSize - 1, refCol:refCol + blockSize - 1);

                currentBlock = double(currentBlock);
                refBlock = double(refBlock);

                % Calculate the Mean Absolute Error 
                MAE = mean(abs(currentBlock(:) - refBlock(:)));

                % Compute the L1 norm of the current motion vector
                currentL1Norm = abs(yOffset) + abs(xOffset);

               % Use Mean of Absolute Error (MAE) as the ME assessment metric. 
               % In case of a tie, choose the block with the smallest motion vector (L1 norm: |𝑥|+|𝑦|, not Euclidean). 
               % If there is still a tie, choose  the  block  with  smallest y.  
               % If there are more than one, choose the one with the smallest x.

                if MAE < minMAE || (MAE == minMAE && (currentL1Norm < bestL1Norm)) || ...
                   (MAE == minMAE && currentL1Norm == bestL1Norm && yOffset < bestVector(1)) || ...
                   (MAE == minMAE && currentL1Norm == bestL1Norm && yOffset == bestVector(1) && xOffset < bestVector(2))

                    minMAE = MAE;
                    bestVector = [yOffset, xOffset];
                    bestL1Norm = currentL1Norm;
                end
            end
        end
    end
end
