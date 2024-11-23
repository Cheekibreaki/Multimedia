% function [approximatedPredictedFrame, predictionModes, vbs_matrix, residualFrame] = vbs_intraPrediction_block_parallel(currentFrame, blockSize, dct_blockSize, baseQP, lambda)
%     % VBS Intra Prediction with Block-Level Wavefront Parallelism
% 
%     [height, width] = size(currentFrame);
%     approximatedPredictedFrame_split = zeros(size(currentFrame), 'double');
%     approximatedPredictedFrame_large = zeros(size(currentFrame), 'double');
%     approximatedPredictedFrame = zeros(size(currentFrame), 'double');
%     approximatedReconstructedFrame_split = zeros(size(currentFrame), 'double');
%     approximatedReconstructedFrame_large = zeros(size(currentFrame), 'double');
%     approximatedReconstructedFrame = zeros(size(currentFrame), 'double');
%     residualFrame_split = zeros(size(currentFrame), 'double');
%     residualFrame_large = zeros(size(currentFrame), 'double');
%     residualFrame = zeros(size(currentFrame), 'double');
%     numBlocksY = ceil(height / blockSize);
%     numBlocksX = ceil(width / blockSize);
%     predictionModes_split = int32(zeros(numBlocksY, numBlocksX));
%     predictionModes_large = int32(zeros(numBlocksY, numBlocksX));
%     predictionModes = int32(zeros(numBlocksY, numBlocksX));
%     vbs_matrix = -1 * ones(numBlocksY, numBlocksX);
% 
%     % Get number of big blocks
%     bigBlockYnum = ceil(numBlocksY / 2);
%     bigBlockXnum = ceil(numBlocksX / 2);
% 
%     % Loop over diagonals for wavefront processing
%     for diagonal = 1:(bigBlockXnum + bigBlockYnum - 1)
%         % Determine the range of big blocks in this diagonal
%         startBigBlockY = max(1, diagonal - bigBlockXnum + 1);
%         endBigBlockY = min(diagonal, bigBlockYnum);
% 
%         % Temp storage for results in this diagonal
%         tempPredictedFrame = zeros(size(currentFrame), 'double');
%         tempReconstructedFrame = zeros(size(currentFrame), 'double');
%         tempResidualFrame = zeros(size(currentFrame), 'double');
%         tempPredictionModes = zeros(numBlocksY, numBlocksX, 'int32');
%         tempVBSMatrix = -1 * ones(numBlocksY, numBlocksX);
% 
%         % Process the current diagonal in parallel
%         for blockIdx = startBigBlockY:endBigBlockY
%             % Local variables for each `parfor` iteration
%             localPredictedFrame = zeros(size(currentFrame), 'double');
%             localReconstructedFrame = zeros(size(currentFrame), 'double');
%             localResidualFrame = zeros(size(currentFrame), 'double');
%             localPredictionModes = zeros(numBlocksY, numBlocksX, 'int32');
%             localVBSMatrix = -1 * ones(numBlocksY, numBlocksX);
% 
%             % Get big block coordinates
%             bigBlockY = blockIdx;
%             bigBlockX = diagonal - bigBlockY + 1;
% 
%             % Translate big block to pixel ranges
%             blockY = 2 * bigBlockY - 1;
%             blockX = 2 * bigBlockX - 1;
%             rowOffset = (blockY - 1) * blockSize + 1;
%             colOffset = (blockX - 1) * blockSize + 1;
%             actualBlockHeight = min(2 * blockSize, height - rowOffset + 1);
%             actualBlockWidth = min(2 * blockSize, width - colOffset + 1);
% 
%             % Extract current block
%             currentBlock = currentFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1);
% 
%             % VBS split estimation
%             [quantized_residualBlock_split, approximatedReconstructed_block_split, approximatedPredictedFrame_split, predictionModes_split, approximatedReconstructedFrame_split] = ...
%                 VBS_split_estimation(currentFrame, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame, ...
%                 blockY, blockX, blockSize, dct_blockSize, baseQP, numBlocksY, numBlocksX);
% 
%             % VBS large estimation
%             [quantized_residualBlock_large, approximatedReconstructed_block_large, approximatedPredictedFrame_large, predictionModes_large, approximatedReconstructedFrame_large] = ...
%                 VBS_large_estimation(currentFrame, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame, ...
%                 blockY, blockX, blockSize, dct_blockSize, baseQP);
% 
%             % Compute SAD for reconstructed_split
%             SAD_split = sum(sum(abs(double(currentBlock) - double(approximatedReconstructed_block_split))));
% 
%             % Compute SAD for reconstructed_large
%             SAD_large = sum(sum(abs(double(currentBlock) - double(approximatedReconstructed_block_large))));
% 
%             % Entropy encode (placeholders)
%             [encodedMDiff_large, encodedResidues_large] = entropyEncode(true, [], predictionModes_large(blockY:blockY+1, blockX:blockX+1), quantized_residualBlock_large);
%             [encodedMDiff_split, encodedResidues_split] = entropyEncode(true, [], predictionModes_split(blockY:blockY+1, blockX:blockX+1), quantized_residualBlock_split);
% 
%             % Calculate rate (R) for large and split blocks
%             total_bits_large = numel(encodedMDiff_large) + numel(encodedResidues_large);
%             total_bits_split = numel(encodedMDiff_split) + numel(encodedResidues_split);
% 
%             R_large = total_bits_large / (total_bits_large + total_bits_split);
%             R_split = total_bits_split / (total_bits_large + total_bits_split);
% 
%             D_large = SAD_large / (SAD_large + SAD_split);
%             D_split = SAD_split / (SAD_large + SAD_split);
% 
%             J_large = D_large + lambda * R_large;
%             J_split = D_split + lambda * R_split;
% 
%             % Select better RD cost
%             if J_large < J_split
%                 localPredictedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = approximatedPredictedFrame_large(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1);
%                 localReconstructedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = approximatedReconstructed_block_large;
%                 localResidualFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = quantized_residualBlock_large;
%                 localPredictionModes(blockY:blockY+1, blockX:blockX+1) = predictionModes_large(blockY:blockY+1, blockX:blockX+1);
%                 localVBSMatrix(blockY:blockY+1, blockX:blockX+1) = 0;
%             else
%                 localPredictedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = approximatedPredictedFrame_split(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1);
%                 localReconstructedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = approximatedReconstructed_block_split;
%                 localResidualFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = quantized_residualBlock_split;
%                 localPredictionModes(blockY:blockY+1, blockX:blockX+1) = predictionModes_split(blockY:blockY+1, blockX:blockX+1);
%                 localVBSMatrix(blockY:blockY+1, blockX:blockX+1) = 1;
%             end
% 
%             % Aggregate into temporary arrays
%             tempPredictedFrame = tempPredictedFrame + localPredictedFrame;
%             tempReconstructedFrame = tempReconstructedFrame + localReconstructedFrame;
%             tempResidualFrame = tempResidualFrame + localResidualFrame;
%             tempPredictionModes = tempPredictionModes + localPredictionModes;
%             tempVBSMatrix = tempVBSMatrix + localVBSMatrix;
%         end
% 
%         % Update shared variables after wavefront
%         approximatedPredictedFrame = approximatedPredictedFrame + tempPredictedFrame;
%         approximatedReconstructedFrame = approximatedReconstructedFrame + tempReconstructedFrame;
%         residualFrame = residualFrame + tempResidualFrame;
%         predictionModes = predictionModes + tempPredictionModes;
%         vbs_matrix = vbs_matrix + tempVBSMatrix;
%     end
% end

function [approximatedPredictedFrame, predictionModes, vbs_matrix, residualFrame] = ...
    vbs_intraPrediction_block_parallel(currentFrame, blockSize, dct_blockSize, baseQP, lambda)
    % VBS Intra Prediction with spmd-based Parallelism

    [height, width] = size(currentFrame);
    approximatedPredictedFrame = zeros(size(currentFrame), 'double');
    approximatedReconstructedFrame = zeros(size(currentFrame), 'double');
    residualFrame = zeros(size(currentFrame), 'double');
    numBlocksY = ceil(height / blockSize);
    numBlocksX = ceil(width / blockSize);
    predictionModes = int32(zeros(numBlocksY, numBlocksX));
    vbs_matrix = -1 * ones(numBlocksY, numBlocksX);

    % Divide work among two workers
    spmd(2)
        localPredictedFrame = zeros(size(currentFrame), 'double');
        localReconstructedFrame = zeros(size(currentFrame), 'double');
        localResidualFrame = zeros(size(currentFrame), 'double');
        localPredictionModes = int32(zeros(numBlocksY, numBlocksX));
        localVBSMatrix = -1 * ones(numBlocksY, numBlocksX);

        % Determine which rows each worker processes
        if spmdIndex == 1
            rowStart = 1;
        else
            rowStart = 2;
        end

        % Process rows assigned to this worker
        for blockY = rowStart:2:numBlocksY
            for blockX = 1:numBlocksX
                % Extract the block for RD cost calculation
                rowOffset = (blockY - 1) * blockSize + 1;
                colOffset = (blockX - 1) * blockSize + 1;
                actualBlockHeight = min(blockSize, height - rowOffset + 1);
                actualBlockWidth = min(blockSize, width - colOffset + 1);
                currentBlock = currentFrame(rowOffset:rowOffset + actualBlockHeight - 1, ...
                                            colOffset:colOffset + actualBlockWidth - 1);

                % Get prediction references (vertical/horizontal borders)
                if blockY > 1
                    topReference = localReconstructedFrame(rowOffset-1, colOffset:colOffset+actualBlockWidth-1);
                else
                    topReference = [];
                end
                if blockX > 1
                    leftReference = localReconstructedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset-1);
                else
                    leftReference = [];
                end

                % VBS Split Estimation
                [quantized_residualBlock_split, reconstructedSplit, predictedSplit, modesSplit] = ...
                    VBS_split_estimation(currentBlock, leftReference, topReference, blockSize, ...
                    dct_blockSize, baseQP, lambda);

                % VBS Large Estimation
                [quantized_residualBlock_large, reconstructedLarge, predictedLarge, modesLarge] = ...
                    VBS_large_estimation(currentBlock, leftReference, topReference, blockSize, ...
                    dct_blockSize, baseQP, lambda);

                % Compute SAD for split and large
                SAD_split = sum(abs(double(currentBlock) - double(reconstructedSplit)), 'all');
                SAD_large = sum(abs(double(currentBlock) - double(reconstructedLarge)), 'all');

                % Entropy encode
                [encodedSplit, encodedLarge] = entropyEncode(modesSplit, quantized_residualBlock_split, ...
                                                             modesLarge, quantized_residualBlock_large);

                % Rate-Distortion Cost Calculation
                R_split = numel(encodedSplit);
                R_large = numel(encodedLarge);
                J_split = SAD_split + lambda * R_split;
                J_large = SAD_large + lambda * R_large;

                % Select lower RD cost
                if J_large < J_split
                    localPredictedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = ...
                        predictedLarge;
                    localReconstructedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = ...
                        reconstructedLarge;
                    localResidualFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = ...
                        quantized_residualBlock_large;
                    localPredictionModes(blockY, blockX) = modesLarge;
                    localVBSMatrix(blockY, blockX) = 0; % Large block
                else
                    localPredictedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = ...
                        predictedSplit;
                    localReconstructedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = ...
                        reconstructedSplit;
                    localResidualFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = ...
                        quantized_residualBlock_split;
                    localPredictionModes(blockY, blockX) = modesSplit;
                    localVBSMatrix(blockY, blockX) = 1; % Split block
                end
            end

            % Send reconstructed row to the other worker
            if blockY < numBlocksY
                spmdSend(localReconstructedFrame(rowOffset:rowOffset+actualBlockHeight-1, :), mod(spmdIndex, 2) + 1);
            end

            % Receive reconstructed row from the other worker
            if blockY > 1
                receivedRow = spmdReceive(mod(spmdIndex, 2) + 1);
                localReconstructedFrame((blockY-2)*blockSize+1:(blockY-1)*blockSize, :) = receivedRow;
            end
        end
    end

    % Combine results from both workers
    approximatedPredictedFrame = approximatedPredictedFrame + localPredictedFrame;
    approximatedReconstructedFrame = approximatedReconstructedFrame + localReconstructedFrame;
    residualFrame = residualFrame + localResidualFrame;
    predictionModes = predictionModes + localPredictionModes;
    vbs_matrix = vbs_matrix + localVBSMatrix;
end



function [quantizedResidualBlock, approximatedReconstructedBlock, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame,approximatedresidualFrame] = ...
    VBS_split_estimation(currentFrame, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame, ...
    blockY, blockX, blockSize, dct_blockSize, baseQP, numBlocksY, numBlocksX)

    isLarge = false;
    [height, width] = size(currentFrame);
    baseRowOffset = (blockY - 1) * blockSize + 1;
    baseColOffset = (blockX - 1) * blockSize + 1;

    % Initialize combined blocks
    combinedQuantizedResidualBlock = zeros(2 * blockSize, 2 * blockSize);
    combinedApproximatedReconstructedBlock = zeros(2 * blockSize, 2 * blockSize);
    

    if baseRowOffset == 1 && baseColOffset == 1
        mode_tl = 'mid-gray';
        mode_tr = 'horizontal';
        mode_bl = 'vertical';
        mode_br = 'compare';
    elseif baseRowOffset == 1
        mode_tl = 'horizontal';
        mode_tr = 'horizontal';
        mode_bl = 'compare';
        mode_br = 'compare';
    elseif baseColOffset == 1
        mode_tl = 'vertical';
        mode_tr = 'compare';
        mode_bl = 'vertical';
        mode_br = 'compare';
    else
        mode_tl = 'compare';
        mode_tr = 'compare';
        mode_bl = 'compare';
        mode_br = 'compare';
    end


    % Offsets within the combined block
    offsets = [0, 0; 0, blockSize; blockSize, 0; blockSize, blockSize];
    modes = {'mid-gray', 'horizontal', 'vertical', 'compare'};

    % Process the four sub-blocks
    for idx = 1:4

        if idx == 1
            mode = mode_tl;
        elseif idx == 2
            mode = mode_tr;
        elseif idx == 3
            mode = mode_bl;
        elseif idx == 4
            mode = mode_br;
        end

        subBlockY = blockY + floor((idx - 1) / 2);
        subBlockX = blockX + mod(idx - 1, 2);
        rowOffsetSub = baseRowOffset + offsets(idx, 1);
        colOffsetSub = baseColOffset + offsets(idx, 2);

        if subBlockY <= numBlocksY && subBlockX <= numBlocksX
            [qResidual, aReconstructedBlock, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame,approximatedresidualFrame] = ...
                processBlock(currentFrame, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame, ...
                subBlockY, subBlockX, rowOffsetSub, colOffsetSub, blockSize, dct_blockSize, baseQP, mode, isLarge);

            % Place the sub-blocks into the combined block
            actualSubBlockHeight = size(aReconstructedBlock, 1);
            actualSubBlockWidth = size(aReconstructedBlock, 2);
            rowRange = offsets(idx, 1) + 1:offsets(idx, 1) + actualSubBlockHeight;
            colRange = offsets(idx, 2) + 1:offsets(idx, 2) + actualSubBlockWidth;

            combinedQuantizedResidualBlock(rowRange, colRange) = qResidual;
            combinedApproximatedReconstructedBlock(rowRange, colRange) = aReconstructedBlock;
        end
    end

    quantizedResidualBlock = combinedQuantizedResidualBlock;
    approximatedReconstructedBlock = combinedApproximatedReconstructedBlock;
end

function [quantizedResidualBlock, approximatedReconstructedBlock, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame,approximatedresidualFrame] = ...
    VBS_large_estimation(currentFrame, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame, ...
    blockY, blockX, blockSize, dct_blockSize, baseQP)

    isLarge = true;

    % % Adjust blockY and blockX for large block positioning
    % adjustedBlockY = (blockY - 1) / 2 + 1;
    % adjustedBlockX = (blockX - 1) / 2 + 1;
    % 
    % % Calculate base offsets for large blocks
    % baseRowOffset = (adjustedBlockY - 1) * blockSize + 1;
    % baseColOffset = (adjustedBlockX - 1) * blockSize + 1;

    baseRowOffset = (blockY - 1) * blockSize + 1;
    baseColOffset = (blockX - 1) * blockSize + 1;


    % Determine prediction mode based on position
    if blockY == 1 && blockX == 1
        mode = 'mid-gray';
    elseif blockY == 1
        mode = 'horizontal';
    elseif blockX == 1
        mode = 'vertical';
    else
        mode = 'compare';
    end

    % Call processBlock with adjusted offsets
    [quantizedResidualBlock, approximatedReconstructedBlock, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame,approximatedresidualFrame] = ...
        processBlock(currentFrame, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame, ...
        blockY, blockX, baseRowOffset, baseColOffset, blockSize*2, dct_blockSize*2, baseQP, mode, isLarge);
end

function [quantizedResidualBlock, approximatedReconstructedBlock, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame,approximatedresidualFrame] = processBlock(...
    currentFrame, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame, ...
    blockY, blockX, rowOffset, colOffset, blockSize, dct_blockSize, baseQP, mode, isLarge)

    [height, width] = size(currentFrame);
    actualBlockHeight = min(blockSize, height - rowOffset + 1);
    actualBlockWidth = min(blockSize, width - colOffset + 1);

    % Ensure block dimensions are valid
    if actualBlockHeight <= 0 || actualBlockWidth <= 0
        % Skip processing this block
        quantizedResidualBlock = [];
        approximatedReconstructedBlock = [];
        return;
    end

    % Initialize prediction
    predBlock = zeros(actualBlockHeight, actualBlockWidth);

    % Prediction logic
    switch mode
        case 'mid-gray'
            predBlock(:) = 128;
            predictionMode = 2;
        case 'horizontal'
            if colOffset > 1
                leftPixels = approximatedReconstructedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset - 1);
                predBlock = repmat(leftPixels, 1, actualBlockWidth);
                predictionMode = 0;
            else
                predBlock(:) = 128;
                predictionMode = 2;
            end
        case 'vertical'
            if rowOffset > 1
                topPixels = approximatedReconstructedFrame(rowOffset - 1, colOffset:colOffset+actualBlockWidth-1);
                predBlock = repmat(topPixels, actualBlockHeight, 1);
                predictionMode = 1;
            else
                predBlock(:) = 128;
                predictionMode = 2;
            end
        case 'compare'
            predictions = {};
            modes = [];
            errors = [];

            % Horizontal Prediction
            if colOffset > 1
                leftPixels = approximatedReconstructedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset - 1);
                horizPred = repmat(leftPixels, 1, actualBlockWidth);
                predictions{end+1} = horizPred;
                modes(end+1) = 0;
                errors(end+1) = mean2(abs(double(currentFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1)) - double(horizPred)));
            end

            % Vertical Prediction
            if rowOffset > 1
                topPixels = approximatedReconstructedFrame(rowOffset - 1, colOffset:colOffset+actualBlockWidth-1);
                vertPred = repmat(topPixels, actualBlockHeight, 1);
                predictions{end+1} = vertPred;
                modes(end+1) = 1;
                errors(end+1) = mean2(abs(double(currentFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1)) - double(vertPred)));
            end

            if ~isempty(errors)
                [~, idx] = min(errors);
                predBlock = predictions{idx};
                predictionMode = modes(idx);
            else
                predBlock(:) = 128;
                predictionMode = 2;
            end
        otherwise
            error('Unknown prediction mode.');
    end

    % Update predicted frame and prediction modes
    approximatedPredictedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = predBlock;
    predictionModes(blockY, blockX) = predictionMode;

    if isLarge
        predictionModes(blockY:blockY+1, blockX:blockX+1) = predictionMode;
    end

    % Compute residual
    residualBlock = currentFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) - predBlock;

    if isempty(residualBlock)
        % Skip processing if residualBlock is empty
        quantizedResidualBlock = [];
        approximatedReconstructedBlock = [];
        return;
    end

    % Quantize and inverse quantize the residual
    quantizedResidualBlock = quantization(residualBlock, dct_blockSize, actualBlockHeight, actualBlockWidth, baseQP);
    approximatedresidualBlock = invquantization(quantizedResidualBlock, dct_blockSize, actualBlockHeight, actualBlockWidth, baseQP);
    
    % Reconstruct the block
    approximatedReconstructedBlock = predBlock + approximatedresidualBlock;
    approximatedReconstructedBlock = double(max(0, min(255, approximatedReconstructedBlock)));
    

    approximatedresidualFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = approximatedresidualBlock;
    % Update the reconstructed frame
    approximatedReconstructedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = approximatedReconstructedBlock;
end
