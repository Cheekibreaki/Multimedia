function [motionVectors, avgMAE, vbs_matrix] = vbs_motionEstimation(currentFrame, originalReferenceFrames, interpolatedReferenceFrames, blockSize, searchRange, dct_blockSize, QP, lambda, FMEEnable, FastME)
    % Get the dimensions of the frame
    [height, width] = size(currentFrame);

    % Initialize the motion vector array and VBS matrix
    numBlocksX = width / blockSize;
    numBlocksY = height / blockSize;
    motionVectors = zeros(numBlocksY, numBlocksX, 3);  % Stores motion vectors for each block.
    vbs_matrix = -1 * ones(numBlocksY, numBlocksX);
    totalMAE = 0;  % To keep track of the total MAE across all blocks

    previous_motion_vector_block = zeros(1, 1, 3);

    for blockY = 1:2:numBlocksY
        for blockX = 1:2:numBlocksX
            if FMEEnable
                referenceFrames = interpolatedReferenceFrames;
            else
                referenceFrames = originalReferenceFrames;
            end

            % Single large block: size is 2 * blockSize (i.e., twice the size in both dimensions)
            currentBlockSize = blockSize * 2;
            rowOffset = (blockY - 1) * blockSize + 1;
            colOffset = (blockX - 1) * blockSize + 1;

            % Extract the current block from the current frame
            currentBlock = currentFrame(rowOffset:rowOffset + currentBlockSize - 1, ...
                                        colOffset:colOffset + currentBlockSize - 1);

            % Compute motion estimation for large block and split blocks
            [motionVector_block_large, total_minMAE_large] = compute_motionVector_block(currentBlock, currentBlockSize, originalReferenceFrames, interpolatedReferenceFrames, rowOffset, colOffset, searchRange, previous_motion_vector_block, true, FMEEnable, FastME);
            [motionVector_block_split, total_minMAE_split] = compute_motionVector_block(currentBlock, currentBlockSize, originalReferenceFrames, interpolatedReferenceFrames, rowOffset, colOffset, searchRange, previous_motion_vector_block, false, FMEEnable, FastME);

            % Compute predicted frames for large and split blocks
            predictedFrame_block_large = compute_predictedFrame_block(referenceFrames, motionVector_block_large, rowOffset, colOffset, blockSize);
            predictedFrame_block_split = compute_predictedFrame_block(referenceFrames, motionVector_block_split, rowOffset, colOffset, blockSize);

            % Compute residuals
            Residuals_block_large = double(currentBlock) - double(predictedFrame_block_large);
            Residuals_block_split = double(currentBlock) - double(predictedFrame_block_split);

            % Construct vbs_matrix for large block (all zeros indicating large blocks)
            vbs_matrix_large = zeros(currentBlockSize / dct_blockSize, currentBlockSize / dct_blockSize);

            % Quantize residuals using vbs_matrix for large block
            quantizedResiduals_large = quantization(Residuals_block_large, dct_blockSize, currentBlockSize, currentBlockSize, QP, vbs_matrix_large);

            % Construct vbs_matrix for split blocks (all ones indicating small blocks)
            vbs_matrix_split = ones(currentBlockSize / dct_blockSize, currentBlockSize / dct_blockSize);

            % Quantize residuals using vbs_matrix for split blocks
            quantizedResiduals_split = quantization(Residuals_block_split, dct_blockSize, currentBlockSize, currentBlockSize, QP, vbs_matrix_split);

            % Inverse quantization
            compresiduals_large = invquantization_block(quantizedResiduals_large, dct_blockSize, currentBlockSize, currentBlockSize, QP, vbs_matrix_large);
            compresiduals_split = invquantization_block(quantizedResiduals_split, dct_blockSize, currentBlockSize, currentBlockSize, QP, vbs_matrix_split);

            % Reconstruct frames
            reconstructed_large = predictedFrame_block_large + compresiduals_large;
            reconstructed_split = predictedFrame_block_split + compresiduals_split;

            % Compute SAD (Sum of Absolute Differences)
            SAD_large = sum(sum(abs(double(currentBlock) - double(reconstructed_large))));
            SAD_split = sum(sum(abs(double(currentBlock) - double(reconstructed_split))));

            % Differential encoding for motion vectors
            [MDiffMV_large, previous_motion_vector_block_large] = diffEncoding_block(motionVector_block_large, 'mv', previous_motion_vector_block);
            [MDiffMV_split, previous_motion_vector_block_split] = diffEncoding_block(motionVector_block_split, 'mv', previous_motion_vector_block);

            % Entropy encoding
            [encodedMDiff_large, ~, encodedResidues_large] = entropyEncode(false, MDiffMV_large, [], quantizedResiduals_large);
            [encodedMDiff_split, ~, encodedResidues_split] = entropyEncode(false, MDiffMV_split, [], quantizedResiduals_split);

            % Calculate rate (R) for large and split blocks
            total_bits_large = numel(encodedMDiff_large) + numel(encodedResidues_large);
            total_bits_split = numel(encodedMDiff_split) + numel(encodedResidues_split);

            total_bits = total_bits_large + total_bits_split;  % Total bits for normalization

            R_large = total_bits_large / total_bits;
            R_split = total_bits_split / total_bits;

            % Normalize SAD for distortion (D)
            total_SAD = SAD_large + SAD_split;
            D_large = SAD_large / total_SAD;
            D_split = SAD_split / total_SAD;

            % Calculate RD cost
            J_large = D_large + lambda * R_large;
            J_split = D_split + lambda * R_split;

            % Choose the motion vector block with the lower RD cost
            if J_large < J_split
                motionVector_block = motionVector_block_large;
                total_minMAE = total_minMAE_large;
                vbs_matrix(blockY:blockY+1, blockX:blockX+1) = 0;
                previous_motion_vector_block = previous_motion_vector_block_large;
            else
                motionVector_block = motionVector_block_split;
                total_minMAE = total_minMAE_split;
                vbs_matrix(blockY:blockY+1, blockX:blockX+1) = 1;
                previous_motion_vector_block = previous_motion_vector_block_split;
            end

            % Store the motion vectors in the corresponding positions
            motionVectors(blockY:blockY+1, blockX:blockX+1, :) = motionVector_block;

            % Update the total MAE
            totalMAE = totalMAE + total_minMAE;
        end
    end

    % Calculate the average MAE across all blocks
    avgMAE = totalMAE / (numBlocksX * numBlocksY / 4);
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




function [motionVector_block,total_minMAE] = compute_motionVector_block(currentBlock, currentBlockSize, originalReferenceFrames, interpolatedReferenceFrames, rowOffset, colOffset, searchRange, predictedMV, isLarge,FMEEnable, FastME)
    
    % Initialize variables for the best match
    bestVector = [0, 0];
    bestRefIdx = 0;
    minMAE = inf;
    bestL1Norm = inf;
    total_minMAE = 0;
    % Split current block into 4 pieces and find the best match for each
    motionVector_block = zeros(2, 2, 3);
    subBlockSize = currentBlockSize / 2;
    predictedMV_new = reshape(predictedMV(:,:,1:2), 1, 2); % Explicitly reshape to 1x2
    predictedMV = predictedMV_new;
    if FMEEnable
        referenceFrames = interpolatedReferenceFrames;
    else
        referenceFrames = originalReferenceFrames;
    end
    if isLarge
        
        % Check all reference frames to find the best match
        for refIdx = 1:length(referenceFrames)
            referenceFrame = referenceFrames{refIdx};
            if FMEEnable
                if FastME
                % If fast ME is enabled
                    [vector, mae, L1Norm] = findBestMatchFastFraction(currentBlock, referenceFrame, rowOffset,  colOffset, currentBlockSize, searchRange, predictedMV);
                else
                % If fast ME is NOT enabled
                    [vector, mae, L1Norm] = findBestMatchFractionalPixel(currentBlock, referenceFrame,rowOffset,  colOffset, currentBlockSize, searchRange);
                end
            else
            % If fractional ME is NOT enabled
                if FastME
                % If fast ME is enabled
                    [vector, mae, L1Norm] = findBestMatchFast(currentBlock, referenceFrame,rowOffset,  colOffset, currentBlockSize,searchRange, predictedMV);
                else
                % If fast ME is NOT enabled
                    [vector, mae, L1Norm] = findBestMatchFullPixel(currentBlock, referenceFrame, rowOffset,  colOffset, currentBlockSize, searchRange);
                end

            end
            
           
    
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
                    if FMEEnable
                        if FastME
                        % If fast ME is enabled
                        
                            [vector, mae, L1Norm] = findBestMatchFastFraction(subBlock, referenceFrame, rowOffset + subBlockY - 1,  colOffset + subBlockX - 1, subBlockSize, searchRange, predictedMV);
                        else
                        % If fast ME is NOT enabled
                            [vector, mae, L1Norm] = findBestMatchFractionalPixel(subBlock, referenceFrame, rowOffset + subBlockY - 1,  colOffset + subBlockX - 1,  subBlockSize, searchRange);
                        end
                    else
                    % If fractional ME is NOT enabled
                        if FastME
                        % If fast ME is enabled
                            [vector, mae, L1Norm] = findBestMatchFast(subBlock, referenceFrame,rowOffset + subBlockY - 1,  colOffset + subBlockX - 1,  subBlockSize,searchRange, predictedMV);
                        else
                        % If fast ME is NOT enabled
                            [vector, mae, L1Norm] = findBestMatchFullPixel(subBlock, referenceFrame,rowOffset + subBlockY - 1,  colOffset + subBlockX - 1, subBlockSize, searchRange);
                        end
        
                    end
                    
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

