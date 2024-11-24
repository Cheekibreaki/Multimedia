function [approximatedPredictedFrame, predictionModes, vbs_matrix, residualFrame] = ...
    vbs_intraPrediction_Mode2(currentFrame, blockSize, dct_blockSize, baseQP, lambda)
    % VBS Intra Prediction with spmd-based Parallelism
    % Each worker processes alternating rows (odd/even). This ensures that
    % before processing a block, the necessary data from previous row
    % (vertical prediction) and from blocks in the same row (horizontal
    % prediction) is available.
    [height, width] = size(currentFrame);
    approximatedPredictedFrame = zeros(size(currentFrame), 'double');
    approximatedReconstructedFrame = zeros(size(currentFrame), 'double');
    residualFrame = zeros(size(currentFrame), 'double');
    numBlocksY = ceil(height / blockSize);
    numBlocksX = ceil(width / blockSize);
    predictionModes = int32(zeros(numBlocksY, numBlocksX));
    vbs_matrix = -1 * ones(numBlocksY, numBlocksX);
     % Get number of big blocks
    bigBlockYnum = ceil(numBlocksY / 2);
    bigBlockXnum = ceil(numBlocksX / 2);

    % Divide work among two workers
    % Worker 1: Process odd rows and send the block to worker 2 as soon as the block is processed. It would start processing the next odd row when receiving the row above it from worker 2.
    % Worker 2: When it gets the block above it from worker 1, it can start processing the even rows.After the whole row is processed, the row is send to worker 1. 
    spmd(2)
        
    % Debugging inside spmd
    % if spmdIndex == 1
    %     disp('Worker 1 - Current vbs_matrix:');
    %     disp(vbs_matrix);
    % elseif spmdIndex == 2
    %     disp('Worker 2 - Current vbs_matrix:');
    %     disp(vbs_matrix);
    % end
       
            localPredictedFrame = zeros(size(currentFrame), 'double');
            localReconstructedFrame = zeros(size(currentFrame), 'double');
            localResidualFrame = zeros(size(currentFrame), 'double');
            localPredictionModes = int32(zeros(numBlocksY, numBlocksX));
            localVBSMatrix = -1 * ones(numBlocksY, numBlocksX);
            
             if spmdIndex == 1
                % Worker 1 processes odd rows 
                rowStart = 1;
            else
                % Worker 2 processes even rows 
                rowStart = 2;
            end
        % Process rows assigned to this worker
            for bigBlockY = rowStart:2:bigBlockYnum
               
                % If not the first row, worker 1 receive reconstructed row from worker 2
                if spmdIndex ==1 
                    if bigBlockY > 1
                    reconstructedRow = spmdReceive;
                    localReconstructedFrame((bigBlockY-2)*2*blockSize+1 : (bigBlockY-1)*2*blockSize, :) = reconstructedRow;
                    end
                end

                for bigBlockX = 1:bigBlockXnum
                    
                    % Worker 2 receives the block directly above it from worker 1
                    if spmdIndex ==2
                        receivedBlock = spmdReceive(1); 
                        localReconstructedFrame((bigBlockY -2)*(2*blockSize)+1 : (bigBlockY -1)*(2*blockSize),...
                                                ((bigBlockX -1)*(2*blockSize)+1 : bigBlockX*2*blockSize)) = receivedBlock;
                    end   
                    
                    % Extract the block for RD cost calculation
                    % Translate from big block idx to small block idx
                    blockY = 2 * bigBlockY -1;
                    blockX = 2 * bigBlockX -1;
    
                    % Translate from small block idx to pixel values
                    rowOffset = (blockY - 1) * blockSize + 1;
                    colOffset = (blockX - 1) * blockSize + 1;
                    actualBlockHeight = 2 * blockSize;
                    actualBlockWidth = 2 * blockSize;
    
                    currentBlock = currentFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1);
    
                    % VBS split estimation
                 [quantized_residualBlock_split, approximatedReconstructed_block_split, approximatedPredictedFrame_split, predictionModes_split, approximatedReconstructedFrame_split] = ...
                     VBS_split_estimation(currentFrame, localPredictedFrame, localPredictionModes, localReconstructedFrame, ...
                     blockY, blockX, blockSize, dct_blockSize, baseQP, numBlocksY, numBlocksX);
                   
                 % VBS large estimation
                [quantized_residualBlock_large, approximatedReconstructed_block_large, approximatedPredictedFrame_large, predictionModes_large, approximatedReconstructedFrame_large] = ...
                    VBS_large_estimation(currentFrame, localPredictedFrame, localPredictionModes, localReconstructedFrame, ...
                    blockY, blockX, blockSize, dct_blockSize, baseQP);
                    
    
                  % Compute SAD for reconstructed_split
                    SAD_split = sum(sum(abs(double(currentBlock) - double(approximatedReconstructed_block_split))));
                  % Compute SAD for reconstructed_large
                    SAD_large = sum(sum(abs(double(currentBlock) - double(approximatedReconstructed_block_large))));
                    
                     [encodedMDiff_large, encodedResidues_large] = entropyEncode(true, [], predictionModes_large(blockY:blockY+1, blockX:blockX+1), quantized_residualBlock_large);
                     [encodedMDiff_split, encodedResidues_split] = entropyEncode(true, [], predictionModes_split(blockY:blockY+1, blockX:blockX+1), quantized_residualBlock_split);
        
                    % Rate-Distortion Cost Calculation
                    % Calculate rate (R) for large and split blocks
                    total_bits_large = numel(encodedMDiff_large) + numel(encodedResidues_large);
                    total_bits_split = numel(encodedMDiff_split) + numel(encodedResidues_split);
        
                    R_large = total_bits_large / (total_bits_large + total_bits_split);
                    R_split = total_bits_split / (total_bits_large + total_bits_split);
        
                    D_large = SAD_large / (SAD_large + SAD_split);
                    D_split = SAD_split / (SAD_large + SAD_split);
        
                    J_large = D_large + lambda * R_large;
                    J_split = D_split + lambda * R_split;
        
        
                    if J_large < J_split
                        localReconstructedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = approximatedReconstructed_block_large;
                        localPredictedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = approximatedPredictedFrame_large(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1);
                        localResidualFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = quantized_residualBlock_large;
                        localPredictionModes(blockY:blockY+1, blockX:blockX+1) = predictionModes_large(blockY:blockY+1, blockX:blockX+1);
                        localVBSMatrix(blockY:blockY+1, blockX:blockX+1) = 0;
                    else
                        localReconstructedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = approximatedReconstructed_block_split;
                        localPredictedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = approximatedPredictedFrame_split(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1);
                        localResidualFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = quantized_residualBlock_split;
                        localPredictionModes(blockY:blockY+1, blockX:blockX+1) = predictionModes_split(blockY:blockY+1, blockX:blockX+1);
                        localVBSMatrix(blockY:blockY+1, blockX:blockX+1) = 1;
                    end
                    
                    % If Worker 1, send reconstructed block to Worker 2
                    if spmdIndex == 1
                        spmdSend(localReconstructedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1), 2);
                    end

                end
                % If Worker 2 sends reconstructed row to Worker 1
                if spmdIndex ==2
                    spmdSend(localReconstructedFrame((bigBlockY-1)*blockSize*2+1 : bigBlockY*blockSize*2,:), 1);
                end
            end
            % Worker 1 receive the last row from worker 2
            if spmdIndex == 1
             spmdReceive;
            end
    end  
            
    % Combine results from both workers
     for bigBlockY = 1:bigBlockYnum
        if mod(bigBlockY, 2) == 1
            % first worker
            approximatedPredictedFrame((bigBlockY-1)* 2*blockSize +1: bigBlockY * 2 * blockSize,:) = localPredictedFrame{1}((bigBlockY-1)* 2 * blockSize +1: bigBlockY * 2 * blockSize,:);
            residualFrame((bigBlockY-1)* 2 * blockSize +1: bigBlockY * 2 * blockSize,:) = localResidualFrame{1}((bigBlockY-1)* 2 * blockSize +1: bigBlockY * 2 * blockSize,:);
            vbs_matrix(2*bigBlockY-1:2*bigBlockY, :) = localVBSMatrix{1}(2*bigBlockY-1:2*bigBlockY, :);
            predictionModes(2*bigBlockY-1:2*bigBlockY, :) = localPredictionModes{1}(2*bigBlockY-1:2*bigBlockY, :);
        else
            % second worker
           approximatedPredictedFrame((bigBlockY-1)* 2 * blockSize +1: bigBlockY * 2 * blockSize,:) = localPredictedFrame{2}((bigBlockY-1)* 2 * blockSize +1: bigBlockY * 2 * blockSize,:);
           residualFrame((bigBlockY-1)* 2 * blockSize +1: bigBlockY * 2 * blockSize,:) = localResidualFrame{2}((bigBlockY-1)* 2 * blockSize +1: bigBlockY * 2 * blockSize,:);
           vbs_matrix(2*bigBlockY-1:2*bigBlockY, :) = localVBSMatrix{2}(2*bigBlockY-1:2*bigBlockY, :);
           predictionModes(2*bigBlockY-1:2*bigBlockY, :) = localPredictionModes{2}(2*bigBlockY-1:2*bigBlockY, :);
        end
     end

    
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
