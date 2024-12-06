function [approximatedPredictedFrame, predictionModes, vbs_matrix,residualFrame,final_encodedResidues] = vbs_intraPrediction(currentFrame, blockSize, dct_blockSize, baseQP,lambda,RCflag,per_block_row_budget, bitCountPerRow)
    addpath('../Utils');  % For utils functions
    [height, width] = size(currentFrame);
    approximatedPredictedFrame_split = zeros(size(currentFrame), 'double');
    approximatedPredictedFrame_large = zeros(size(currentFrame), 'double');
    approximatedPredictedFrame = zeros(size(currentFrame), 'double');
    approximatedReconstructedFrame_split = zeros(size(currentFrame), 'double');
    approximatedReconstructedFrame_large = zeros(size(currentFrame), 'double');
    approximatedReconstructedFrame = zeros(size(currentFrame), 'double');
    numBlocksY = ceil(height / blockSize);
    numBlocksX = ceil(width / blockSize);
    residualFrame = zeros(size(currentFrame), 'double');
    predictionModes_split = int32(zeros(numBlocksY, numBlocksX));
    predictionModes_large = int32(zeros(numBlocksY, numBlocksX));
    predictionModes = int32(zeros(numBlocksY, numBlocksX));
    vbs_matrix = -1 * ones(numBlocksY, numBlocksX);

    final_encodedResidues = []
    row_bits_used = per_block_row_budget;
    for blockY = 1:2:numBlocksY
         if RCflag
            next_row_budget = per_block_row_budget + (per_block_row_budget - row_bits_used);
            baseQP = findCorrectQP(next_row_budget,bitCountPerRow);
            baseQP = 0
            row_bits_used = 0;
        end
        for blockX = 1:2:numBlocksX
            % Extract the current block for RD cost calculation
            rowOffset = (blockY - 1) * blockSize + 1;
            colOffset = (blockX - 1) * blockSize + 1;
            actualBlockHeight = min(2 * blockSize, height - rowOffset + 1);
            actualBlockWidth = min(2 * blockSize, width - colOffset + 1);
            currentBlock = currentFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1);

            % VBS split estimation
            [quantized_residualBlock_split, approximatedReconstructed_block_split, approximatedPredictedFrame_split, predictionModes_split, approximatedReconstructedFrame_split] = ...
                VBS_split_estimation(currentFrame, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame, ...
                blockY, blockX, blockSize, dct_blockSize, baseQP, numBlocksY, numBlocksX);

            % VBS large estimation
            [quantized_residualBlock_large, approximatedReconstructed_block_large, approximatedPredictedFrame_large, predictionModes_large, approximatedReconstructedFrame_large] = ...
                VBS_large_estimation(currentFrame, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame, ...
                blockY, blockX, blockSize, dct_blockSize, baseQP);
            
            % Compute SAD for reconstructed_split
            SAD_split = sum(sum(abs(double(currentBlock) - double(approximatedReconstructed_block_split))));

            % Compute SAD for reconstructed_large
            SAD_large = sum(sum(abs(double(currentBlock) - double(approximatedReconstructed_block_large))));

            % Placeholder for entropyEncode function (you need to implement this)
            % Assuming entropyEncode returns encoded data and its length
            non_important1 = [];
            [nonimportant1, encodedMDiff_large, encodedResidues_large] = entropyEncode(true, non_important1,predictionModes_large(blockY:blockY+1, blockX:blockX+1), quantized_residualBlock_large);
            [nonimportant1,encodedMDiff_split, encodedResidues_split] = entropyEncode(true, non_important1,predictionModes_split(blockY:blockY+1, blockX:blockX+1), quantized_residualBlock_split);



            % Calculate rate (R) for large and split blocks
            total_bits_large = numel(encodedMDiff_large) + numel(encodedResidues_large);
            total_bits_split = numel(encodedMDiff_split) + numel(encodedResidues_split);

            R_large = total_bits_large/total_bits_large+total_bits_split;
            R_split = total_bits_split/total_bits_large+total_bits_split;
        
            D_large = SAD_large/SAD_split+SAD_large;
            D_split = SAD_split/SAD_split+SAD_large;


            % Calculate RD cost

            J_large = D_large + lambda * R_large;
            J_split = D_split + lambda * R_split;

            % Choose the block with the lower RD cost
            if J_large < J_split
                approximatedReconstructedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = ...
                    approximatedReconstructed_block_large;
                
                approximatedPredictedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = ...
                    approximatedPredictedFrame_large(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1);
                predictionModes(blockY:blockY+1, blockX:blockX+1) = predictionModes_large(blockY:blockY+1, blockX:blockX+1);
                vbs_matrix(blockY:blockY+1, blockX:blockX+1) = 0;
                residualFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = quantized_residualBlock_large;
                encodedResidues = encodedResidues_large;
            else
                approximatedReconstructedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = ...
                    approximatedReconstructed_block_split;
                approximatedPredictedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = ...
                    approximatedPredictedFrame_split(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1);
                predictionModes(blockY:blockY+1, blockX:blockX+1) = predictionModes_split(blockY:blockY+1, blockX:blockX+1);
                vbs_matrix(blockY:blockY+1, blockX:blockX+1) = 1;
                residualFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = quantized_residualBlock_split;
                encodedResidues = encodedResidues_split;
            end
            if RCflag
                encodedResidues_length = length(encodedResidues);
             % Flatten the quantized block using zigzag scan
                row_bits_used = row_bits_used + encodedResidues_length;
           
            end
            final_encodedResidues = [final_encodedResidues, -1 ,baseQP,encodedResidues];
        end
         
    end
    % reconstructedImage = uint8(approximatedPredictedFrame_large);
    %         imshow(reconstructedImage);
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
                subBlockY, subBlockX, rowOffsetSub, colOffsetSub, blockSize, dct_blockSize, baseQP -1, mode, isLarge);

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
