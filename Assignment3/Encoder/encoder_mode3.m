function encoder_mode3(referenceFile, paddedOutputFile, numFrames, width, height, blockSize, searchRange, dct_blockSize, QP, I_Period, nRefFrames,lambda,VBSEnable,FMEEnable,FastME,mode)
% The encoder for mode 3 frame level parallelism is similar to block level
% implementation. Worker 1 processes odd num frames and worker 2 processes
% even num frames. After worker one processes 2 row of blocks, it would
% send the reference blocks to worker 2 immediately. Worker 2 would only
% the reference frame to worker 1 after it processes the whole frame.
    lambda = lambda;
    mode = mode;
    % Parameters for decoder
    params.width = width;                
    params.height = height;                
    params.numFrames = numFrames;              
    params.blockSize = blockSize;                        
    params.dct_blockSize = dct_blockSize;           
    params.QP = QP;                   
    params.nRefFrames = nRefFrames;
    params.VBSEnable = VBSEnable;
    params.FMEEnable = FMEEnable;
    params.FastME = FastME;
    
    % Save the parameters to a MAT-file
    save('../Outputs/headerfile.mat', 'params');

    % Write parameters needed for decoder to header file
    headerFile = fopen('../Outputs/headerfile.mat', 'w');
    fwrite(headerFile, [width, height, numFrames, blockSize, dct_blockSize, QP, nRefFrames,VBSEnable, FMEEnable, FastME], 'int32');
    fclose(headerFile);

    % Open the padded Y only file
    fid = fopen(paddedOutputFile, 'r');

    if fid == -1
        error('Failed to open input file: %s', paddedOutputFile);
    end

    % Initialize storage for frames
    frames = cell(1, numFrames);
    for frameIdx = 1:numFrames
        frames{frameIdx} = fread(fid, [width, height], 'uint8')';
    end
    fclose(fid);

    % if FMEEnable
    %     referenceFrame = 128 * ones(height, width, 'uint8');
    % else
    %     referenceFrame = 128 * ones(2*height - 1, 2*width - 1, 'uint8');
    % end

    % Initialize a cell array to collect data from workers
    frameData = cell(1, numFrames);

    % Frame-level parallelism
    spmd(2)
        % Each worker processes alternating frames
        if spmdIndex == 1
            % Worker 1 processes odd-numbered frames
            frameStart = 1;
        else
            % Worker 2 processes even-numbered frames
            frameStart = 2;
        end

        for frameIdx = frameStart:2:numFrames
            
            if spmdIndex ==1 
                if frameIdx == 1
                    if FMEEnable
                        referenceFrame = 128 * ones(height, width, 'uint8');
                        referenceFrames = {referenceFrame};
                    else 
                        referenceFrame = 128 * ones(2*height - 1, 2*width - 1, 'uint8');
                        referenceFrames = {referenceFrame};
                    end
                else
                    fprintf('Worker %d waiting to receive reference frame for frame %d from Worker %d\n', spmdIndex, frameIdx, 2);
                    referenceFrame = spmdReceive(2);
                    if FMEEnable
                        referenceFrame = interpolateFrame(referenceFrame);
                    end
                    fprintf('Worker %d received reference frame for frame %d from Worker %d\n', spmdIndex, frameIdx, 2);
                    referenceFrames = {referenceFrame};
                    
                end
            end

            if spmdIndex ==2
                referenceFrame = spmdReceive(1);
                fprintf('Frame start 1:Worker %d received reference frame for frame %d from Worker %d\n', spmdIndex, frameIdx, 1);
                referenceFrames = {referenceFrame};
            end


            currentFrame = frames{frameIdx};
    
            %determine the frame type
            if frameIdx == 1
                isIFrame = true;
            elseif frameIdx < I_Period
                isIFrame = false;
            else
                isIFrame = (mod(frameIdx - 1, I_Period) == 0);
            end
             %isIFrame = false;
    
            if isIFrame

               if VBSEnable
                        
                        approximatedPredictedFrame = zeros(size(currentFrame), 'double');
                        
                        approximatedReconstructedFrame = zeros(size(currentFrame), 'double');
                        numBlocksY = ceil(height / blockSize);
                        numBlocksX = ceil(width / blockSize);
                        residualFrame = zeros(size(currentFrame), 'double');
                        currPredictionModes = int32(zeros(numBlocksY, numBlocksX));
                        vbs_matrix = -1 * ones(numBlocksY, numBlocksX);
                    
                        for blockY = 1:2:numBlocksY
                             %dummy receive to prevent blocking worker 1
                           if spmdIndex == 2
                                if blockY ~= numBlocksY-1
                                    fprintf('Worker %d waiting to receive reference frame for frame %d from Worker %d\n', spmdIndex, frameIdx, 1);
                                    referenceFrame = spmdReceive(1);
                                    fprintf('Worker %d received reference frame for frame %d from Worker %d\n', spmdIndex, frameIdx, 1);
                                end
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
                                    VBS_split_estimation(currentFrame, approximatedPredictedFrame, currPredictionModes, approximatedReconstructedFrame, ...
                                    blockY, blockX, blockSize, dct_blockSize, QP, numBlocksY, numBlocksX);
                    
                                % VBS large estimation
                                [quantized_residualBlock_large, approximatedReconstructed_block_large, approximatedPredictedFrame_large, predictionModes_large, approximatedReconstructedFrame_large] = ...
                                    VBS_large_estimation(currentFrame, approximatedPredictedFrame, currPredictionModes, approximatedReconstructedFrame, ...
                                    blockY, blockX, blockSize, dct_blockSize, QP);
                                
                                % Compute SAD for reconstructed_split
                                SAD_split = sum(sum(abs(double(currentBlock) - double(approximatedReconstructed_block_split))));
                    
                                % Compute SAD for reconstructed_large
                                SAD_large = sum(sum(abs(double(currentBlock) - double(approximatedReconstructed_block_large))));
                    
                                % Placeholder for entropyEncode function (you need to implement this)
                                % Assuming entropyEncode returns encoded data and its length
                                non_important1 = [];
                                [encodedMDiff_large, encodedResidues_large] = entropyEncode(true, non_important1,predictionModes_large(blockY:blockY+1, blockX:blockX+1), quantized_residualBlock_large);
                                [encodedMDiff_split, encodedResidues_split] = entropyEncode(true, non_important1,predictionModes_split(blockY:blockY+1, blockX:blockX+1), quantized_residualBlock_split);
                    
                                % Calculate rate (R) for large and split blocks
                                total_bits_large = numel(encodedMDiff_large) + numel(encodedResidues_large);
                                total_bits_split = numel(encodedMDiff_split) + numel(encodedResidues_split);
                    
                                R_large = total_bits_large/(total_bits_large+total_bits_split);
                                R_split = total_bits_split/(total_bits_large+total_bits_split);
                            
                                D_large = SAD_large/(SAD_split+SAD_large);
                                D_split = SAD_split/(SAD_split+SAD_large);
                    
                    
                                % Calculate RD cost
                    
                                J_large = D_large + lambda * R_large;
                                J_split = D_split + lambda * R_split;
                    
                                % Choose the block with the lower RD cost
                                if J_large < J_split
                                    approximatedReconstructedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = ...
                                        approximatedReconstructed_block_large;
                                    
                                    approximatedPredictedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = ...
                                        approximatedPredictedFrame_large(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1);
                                    currPredictionModes(blockY:blockY+1, blockX:blockX+1) = predictionModes_large(blockY:blockY+1, blockX:blockX+1);
                                    vbs_matrix(blockY:blockY+1, blockX:blockX+1) = 0;
                                    residualFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = quantized_residualBlock_large;
                                else
                                    approximatedReconstructedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = ...
                                        approximatedReconstructed_block_split;
                                    approximatedPredictedFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = ...
                                        approximatedPredictedFrame_split(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1);
                                    currPredictionModes(blockY:blockY+1, blockX:blockX+1) = predictionModes_split(blockY:blockY+1, blockX:blockX+1);
                                    vbs_matrix(blockY:blockY+1, blockX:blockX+1) = 1;
                                    residualFrame(rowOffset:rowOffset + actualBlockHeight - 1, colOffset:colOffset + actualBlockWidth - 1) = quantized_residualBlock_split;
                                end
                    
                               
                            end
                             
                            if spmdIndex == 1
                                if frameIdx ~= numFrames
                                    fprintf('Worker %d sending frame  %d to Worker %d\n', spmdIndex, frameIdx, 2);
                                    spmdSend(approximatedReconstructedFrame,2)
                                end
                            end
                        end

                     
                   if spmdIndex == 2
                       if frameIdx ~= numFrames
                           fprintf('Worker %d sending frame  %d to Worker %d\n', spmdIndex, frameIdx, 1);
                           spmdSend(approximatedReconstructedFrame, 1);
                       end    
                   end

                   % if spmdIndex == 2
                   %     msg = spmdReceive(1);
                   % end

                   % Since worker 1 will finish scanning frame early than worker 2,
                   % continue to send to avoid deadlock
                   % if spmdIndex == 1
                   %     if frameIdx ~= numFrames
                   %          spmdSend(approximatedReconstructedFrame, 2);
                   %          spmdSend(approximatedReconstructedFrame, 2);
                   %     end
                   % end


                   MDiffModes = currPredictionModes;
                   reconstructedFrame = approximatedReconstructedFrame;
                   predictedFrame = approximatedPredictedFrame;
                   
               else
               
                    approximatedPredictedFrame = zeros(size(currentFrame), 'double');
                    currPredictionModes = int32(zeros(ceil(height/blockSize), ceil(width/blockSize)));
                    approximatedReconstructedFrame(1:height,1:width) = zeros(size(currentFrame), 'double');
                    
                    
                    for y = 1:blockSize:height
                        %dummy receive to prevent blocking worker 1
                       if spmdIndex == 2
                           
                            if y ~= height - blockSize + 1  % Only receive if not at the last row
                                fprintf('Worker %d waiting to receive reference frame for frame %d from Worker %d\n', spmdIndex, frameIdx, 1);
                                referenceFrame = spmdReceive(1);
                                fprintf('Worker %d received reference frame for frame %d from Worker %d\n', spmdIndex, frameIdx, 1);
                            end
                       end
                        for x = 1:blockSize:width
                            actualBlockHeight = min(blockSize, height-y+1);
                            actualBlockWidth = min(blockSize, width-x+1);
                          
                            if y == 1 && x == 1
                                % For the first block, predict with mid-gray
                                approximatedPredictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = 128;
                                currPredictionModes(1, 1) = 2; % 2 for mid-gray prediction
                            else
                                % Horizontal prediction
                                if x > 1
                                    horizPred = repmat(approximatedReconstructedFrame(y:y+actualBlockHeight-1, x-1), 1, actualBlockWidth);
                                else
                                    horizPred = repmat(128, actualBlockHeight, actualBlockWidth);
                                end
                
                                % Vertical prediction
                                if y > 1
                                    vertPred = repmat(approximatedReconstructedFrame(y-1, x:x+actualBlockWidth-1), actualBlockHeight, 1);
                                else
                                    vertPred = repmat(128, actualBlockHeight, actualBlockWidth);
                                end
                
                                % Calculate MAE for both predictions
                                currentBlock = currentFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1);
                                maeHoriz = mean(abs(double(currentBlock) - double(horizPred)), 'all');
                                maeVert = mean(abs(double(currentBlock) - double(vertPred)), 'all');
                                if maeHoriz <= maeVert
                                    approximatedPredictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = horizPred;
                                    currPredictionModes(ceil(y/blockSize), ceil(x/blockSize)) = 0; % 0 for horizontal
                                else
                                    approximatedPredictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = vertPred;
                                    currPredictionModes(ceil(y/blockSize), ceil(x/blockSize)) = 1; % 1 for vertical
                                end
                            end
                            
                            residualBlock = currentFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) - approximatedPredictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1);
                            quantizedBlock = quantization(residualBlock, dct_blockSize,blockSize,blockSize,QP);
                            approximatedresidualBlock = invquantization(quantizedBlock,dct_blockSize,blockSize,blockSize, ...
                                QP);
                            approximatedReconstructedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = double(max(0, min (255, approximatedresidualBlock + approximatedPredictedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1))));      
                        end
                        if spmdIndex == 1
                            if frameIdx ~= numFrames
                                fprintf('Worker %d sending frame  %d to Worker %d\n', spmdIndex, frameIdx, 2);
                                spmdSend(approximatedReconstructedFrame,2)
                            end
                        end
                    end
                   
                   if spmdIndex == 2
                       if frameIdx ~= numFrames
                           fprintf('Worker %d sending frame  %d to Worker %d\n', spmdIndex, frameIdx, 1);
                           spmdSend(approximatedReconstructedFrame, 1);
                       end    
                   end

                   % if spmdIndex == 2
                   %     msg = spmdReceive(1);
                   % end

                   % Since worker 1 will finish scanning frame early than worker 2,
                   % continue to send to avoid deadlock
                   % if spmdIndex == 1
                   %     if frameIdx ~= numFrames
                   %          spmdSend(approximatedReconstructedFrame, 2);
                   %          spmdSend(approximatedReconstructedFrame, 2);
                   %     end
                   % end

                   reconstructedFrame = approximatedReconstructedFrame;
                   predictedFrame = approximatedPredictedFrame;
                   MDiffModes = diffEncoding(currPredictionModes,'modes');
                   Residuals = double(currentFrame) - double(predictedFrame);
               end
               
              
               if VBSEnable
                   quantizedResiduals = residualFrame;
                    [nonimporatant1,encodedMDiff,encodedResiduals] = entropyEncode(isIFrame, [], MDiffModes, quantizedResiduals, vbs_matrix);
               else
                   quantizedResiduals = quantization(Residuals, dct_blockSize,width,height,QP); 
                   [nonimporatant1,encodedMDiff,encodedResiduals] = entropyEncode(isIFrame, [], MDiffModes, quantizedResiduals);
               end
    
          

               workerData{(frameIdx - spmdIndex + 2) / 2} = struct(...
                    'frameIdx', frameIdx, ...
                    'encodedMDiff', encodedMDiff, ...
                    'encodedResiduals', encodedResiduals, ...
                    'reconstructedFrame',reconstructedFrame,...
                    'predictedFrame',predictedFrame);
              
            else
    
                % Motion estimation
                if VBSEnable
                    % Initialize the motion vector array (for storing motion vectors for each block)
                    numBlocksX = width / blockSize;
                    numBlocksY = height / blockSize;
                    motionVectors = zeros(numBlocksY, numBlocksX, 3);  % Stores motion vector for each block.
                    vbs_matrix = -1 * ones(numBlocksY, numBlocksX);
                    totalMAE = 0;  % To keep track of the total MAE across all blocks
                    approximatedReconstructedFrame(1:height,1:width) = zeros(size(currentFrame), 'double');
                    approximatedPredictedFrame = zeros(size(currentFrame), 'double');
                    
                    if spmdIndex == 1
                        referenceFrames = {referenceFrame};
                    end

                    for blockY = 1:2:numBlocksY
                        previous_motion_vector_block = zeros(1, 1, 3);
                         if spmdIndex == 2
                            if blockY ~= numBlocksY -1
                                fprintf('Worker %d waiting to receive reference frame for frame %d from Worker %d\n', spmdIndex, frameIdx, 1);
                                referenceFrame = spmdReceive(1);
                                fprintf('Worker %d received reference frame for frame %d from Worker %d\n', spmdIndex, frameIdx, 1);
                                if FMEEnable
                                    referenceFrame = interpolateFrame(referenceFrame);
                                end
                                referenceFrames = {referenceFrame};
                            end
                         end

                        for blockX = 1:2:numBlocksX
                            % Single large block: size is 2 * blockSize (i.e., twice the size in both dimensions)
                            currentBlockSize = blockSize * 2;
                            rowOffset = (blockY - 1) * blockSize + 1;
                            colOffset = (blockX - 1) * blockSize + 1;
                
                            % Extract the current block from the current frame
                            currentBlock = currentFrame(rowOffset:rowOffset + currentBlockSize - 1, ...
                                                        colOffset:colOffset + currentBlockSize - 1);
                
                            % Compute motion estimation for large block and split blocks
                            
                            [motionVector_block_large, total_minMAE_large] = compute_motionVector_block(currentBlock, currentBlockSize, referenceFrames, rowOffset, colOffset, searchRange, previous_motion_vector_block, true, FMEEnable, FastME);
                            [motionVector_block_split, total_minMAE_split] = compute_motionVector_block(currentBlock, currentBlockSize, referenceFrames, rowOffset, colOffset, searchRange, previous_motion_vector_block,  false, FMEEnable, FastME);

                            % Compute predicted frames for large and split blocks
                            predictedFrame_block_large = compute_predictedFrame_block(referenceFrames, motionVector_block_large, rowOffset, colOffset, blockSize, FMEEnable);
                            predictedFrame_block_split = compute_predictedFrame_block(referenceFrames, motionVector_block_split, rowOffset, colOffset, blockSize, FMEEnable);
                            
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
                                approximatedReconstructedFrame(rowOffset:rowOffset + currentBlockSize - 1, colOffset:colOffset + currentBlockSize - 1) = double(max(0,min(255,reconstructed_large)));
                                approximatedPredictedFrame(rowOffset:rowOffset + currentBlockSize - 1, colOffset:colOffset + currentBlockSize - 1) = double(max(0,min(255,predictedFrame_block_large)));

                            else
                                motionVector_block = motionVector_block_split;
                                total_minMAE = total_minMAE_split;
                                vbs_matrix(blockY:blockY+1, blockX:blockX+1) = 1;
                                previous_motion_vector_block = previous_motion_vector_block_split;
                                approximatedReconstructedFrame(rowOffset:rowOffset + currentBlockSize - 1, colOffset:colOffset + currentBlockSize - 1) = double(max(0,min(255,reconstructed_split)));
                                approximatedPredictedFrame(rowOffset:rowOffset + currentBlockSize - 1, colOffset:colOffset + currentBlockSize- 1) = double(max(0,min(255,predictedFrame_block_split)));
                            end

                            % Store the motion vectors in the corresponding positions
                            motionVectors(blockY:blockY+1, blockX:blockX+1, :) = motionVector_block;

                            totalMAE = totalMAE + total_minMAE;
                        end

                        if spmdIndex == 1
                            if frameIdx ~= numFrames
                                fprintf('Worker %d sending frame  %d to Worker %d\n', spmdIndex, frameIdx, 2);
                                spmdSend(approximatedReconstructedFrame,2)
                            end
                        end

                    end

                    if spmdIndex == 2
                       if frameIdx ~= numFrames
                           fprintf('Worker %d sending frame  %d to Worker %d\n', spmdIndex, frameIdx, 1);
                           spmdSend(approximatedReconstructedFrame, 1);
                       end    
                    end

                    % Calculate the average MAE across all blocks
                    avgMAE = totalMAE / (numBlocksX * numBlocksY /4);
                    
                    MDiffMV = motionVectors;
                    
                    predictedFrame = approximatedPredictedFrame;
                    reconstructedFrame = approximatedReconstructedFrame;
                    Residuals = double(currentFrame) - double(predictedFrame);   
                    
                else
                    fprintf('Worker %d processing a P frame %d \n', spmdIndex, frameIdx);
                    if spmdIndex == 1
                        referenceFrames = {referenceFrame};
                    end
                    approximatedReconstructedFrame(1:height,1:width) = zeros(size(currentFrame), 'double');
                    approximatedPredictedFrame = zeros(size(currentFrame), 'double');
                    % Initialize the motion vector array (for storing motion vectors for each block)
                    numBlocksX = width / blockSize; 
                    numBlocksY = height / blockSize; 
                    currMotionVectors = zeros(numBlocksY, numBlocksX, 3);  % Stores motion vector for each block.
                
                    totalMAE = 0;  % To keep track of the total MAE across all blocks
                    
                    % Process blocks in raster order (left to right, top to bottom)
                    for row = 1:blockSize:height
                        predictedMV = [0,0];
                        if spmdIndex == 2
                            if row ~= height - blockSize + 1  % Only receive if not at the last row
                                fprintf('Worker %d waiting for reference frame for %d\n', spmdIndex, frameIdx);
                                referenceFrame = spmdReceive(1);
                                fprintf('Worker %d received reference frame for %d\n', spmdIndex, frameIdx);
                                if FMEEnable
                                    referenceFrame = interpolateFrame(referenceFrame);
                                end
                                referenceFrames = {referenceFrame};
                            end
                         end
                        for col = 1:blockSize:width
                            % Extract the current block from the current frame
                            currentBlock = currentFrame(row:row+blockSize-1, col:col+blockSize-1);
                            % Initialize variables for the best match
                            bestRefIdx = 0;
                            minMAE = inf;
                            bestL1Norm = inf;
                                % Check all reference frames to find the best match
                            for refIdx = 1:length(referenceFrames)
                                referenceFrame = referenceFrames{refIdx}; 
                                    % If fractional ME is enabled
                                    if FMEEnable
                                        if FastME
                                        % If fast ME is enabled
                                            [vector, mae, L1Norm] = findBestMatchFastFraction(currentBlock, referenceFrame, row, col, blockSize, searchRange, predictedMV);
                                        else
                                        % If fast ME is NOT enabled
                                            [vector, mae, L1Norm] = findBestMatchFractionalPixel(currentBlock, referenceFrame, row, col, blockSize, searchRange);
                                        end
                                    else
                                    % If fractional ME is NOT enabled
                                        if FastME
                                        % If fast ME is enabled
                                            [vector, mae, L1Norm] = findBestMatchFast(currentBlock, referenceFrame, row, col, blockSize, searchRange, predictedMV);
                                        else
                                        % If fast ME is NOT enabled
                                            [vector, mae, L1Norm] = findBestMatchFullPixel(currentBlock, referenceFrame, row, col, blockSize, searchRange);
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
                              
                
                            dy = bestVector(1);
                            dx = bestVector(2);
                            predictedMV = bestVector;
                            % Convert row and col to block indices
                            blockY = (row - 1) / blockSize + 1;  
                            blockX = (col - 1) / blockSize + 1;  
                
                            % Store the motion vector for this block
                            currMotionVectors(blockY, blockX, 1) = dy;
                            currMotionVectors(blockY, blockX, 2) = dx;
                            currMotionVectors(blockY, blockX, 3) = bestRefIdx - 1;
                            % Add the MAE of this block to the total MAE
                            
                            % Compute the coordinates of the reference block based on the motion vector
                            refRowStart = row + dy;
                            refColStart = col + dx;
                
                            % Extract the reference block from the reference frame
                            if FMEEnable
                                refRowStart = 2*row-1 + dy;
                                refColStart = 2*col-1 + dx;
                                refBlock = referenceFrame(refRowStart:2:(refRowStart + 2* blockSize - 2), refColStart:2:(refColStart + 2 * blockSize - 2));
                            else
                                refBlock = referenceFrame(refRowStart:(refRowStart + blockSize - 1), refColStart:(refColStart + blockSize - 1));
                            end
                            
                            
                            % Place the reference block into the predicted frame
                            approximatedPredictedFrame(row:(row + blockSize - 1), col:(col + blockSize - 1)) = refBlock;
                            %double(max(0,min(255,refBlock)));
                            

                            % Compute residuals
                            residuals_block = double(currentBlock) - double(refBlock);
                            quantized_block = quantization(residuals_block, dct_blockSize, blockSize, blockSize, QP);
                            approximated_residualBlock = invquantization(quantized_block,dct_blockSize,blockSize,blockSize,QP);
                            approximatedReconstructedBlock = approximated_residualBlock + refBlock;
                            approximatedReconstructedFrame(row:(row + blockSize - 1), col:(col + blockSize - 1)) = approximatedReconstructedBlock;
                            %double(max(0,min(255,approximatedReconstructedBlock)));
             
                            totalMAE = totalMAE + minMAE;

                            
                        end

                        if spmdIndex == 1
                            if frameIdx ~= numFrames
                                spmdSend(approximatedReconstructedFrame,2)
                            end
                        end
                    end
                
                    % Calculate the average MAE across all blocks
                    avgMAE = totalMAE / (numBlocksX * numBlocksY);
                    % currMotionVectors
                    MDiffMV = diffEncoding(currMotionVectors,'mv');
                    reconstructedFrame = approximatedReconstructedFrame;

                    if spmdIndex == 2
                       if frameIdx ~= numFrames
                           fprintf('Worker %d sending frame %d to worker 1\n', spmdIndex, frameIdx);
                           spmdSend(approximatedReconstructedFrame, 1);
                       end    
                    end
          
                    predictedFrame = approximatedPredictedFrame;
                    Residuals = double(currentFrame) - double(predictedFrame);    
                   
                    % if frameIdx == 8
                    % 
                    %     fprintf('Residual value for frame 8',Residuals)
                    % end

                    % Check if the matrix contains -1 or 0 
                    if any(Residuals(:) == -1 | Residuals(:) == 0) 
                        disp(Residuals); 
                    end
                end

                if VBSEnable
                     quantizedResiduals = quantization(Residuals, dct_blockSize,width,height,QP,vbs_matrix); 
                    [encodedMDiff,nonimporatant1,encodedResiduals] = entropyEncode(isIFrame, MDiffMV, [], quantizedResiduals,vbs_matrix);
                else
                    quantizedResiduals = quantization(Residuals, dct_blockSize,width,height,QP); 
                    [encodedMDiff,nonimporatant1,encodedResiduals] = entropyEncode(isIFrame, MDiffMV, [], quantizedResiduals);
                end
    
               % Store worker data for this frame
                workerData{(frameIdx - spmdIndex + 2) / 2} = struct(...
                    'frameIdx', frameIdx, ...
                    'encodedMDiff', encodedMDiff, ...
                    'encodedResiduals', encodedResiduals, ...
                    'reconstructedFrame',reconstructedFrame,...
                    'predictedFrame',predictedFrame);

            end
        end
    end

    % Combine data from all workers
    allWorkerData = [workerData{1}, workerData{2}];
    
    
    % Extract frame indices
    frameIndices = cellfun(@(data) data.frameIdx, allWorkerData);
    
    % Sort the data by frameIdx
    [~, sortedOrder] = sort(frameIndices);
    sortedWorkerData = allWorkerData(sortedOrder);
    reconstructedYUVFile = fopen(referenceFile, 'w');
    if reconstructedYUVFile == -1
        error('Failed to open YUV file for writing reconstructed frames.');
    end
    
    predictedFile = '../Outputs/predictedFrames.yuv';
    predictedYUVFile = fopen(predictedFile,'w');
    % Save data in the required format
    for i = 1:length(sortedWorkerData)
        frameData = sortedWorkerData{i};
        frameIdx = frameData.frameIdx;
    
        % Save motion vector file
        motionVectorFile = sprintf('../Outputs/MDiff_frame_%d.mat', frameIdx);
        encodedMDiff = frameData.encodedMDiff;
        save(motionVectorFile, 'encodedMDiff');
    
        % Save residuals file
        residualFile = sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx);
        encodedResiduals = frameData.encodedResiduals;
        save(residualFile, 'encodedResiduals');

        % Save reconstructed frame to YUV file
        reconstructedFrame = frameData.reconstructedFrame;
        fwrite(reconstructedYUVFile, reconstructedFrame', 'uint8'); % Write in YUV format

        predictedFrame = frameData.predictedFrame;
        fwrite(predictedYUVFile, predictedFrame','uint8');
    end

    fclose(reconstructedYUVFile);
    fclose(predictedYUVFile);
       
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



function predictedBlock = compute_predictedFrame_block(referenceFrames, motionVector_block, rowOffset, colOffset, blockSize, FMEEnable)
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
            [height, width] = size(referenceFrame); 
            % fprintf('size of referenceFrames %d,%d\n',height,width);
            % Calculate the position of the reference sub-block based on the motion vector
            refRowStart = subBlockRowOffset + mvY;
            refColStart = subBlockColOffset + mvX;
            
            % % Ensure the reference sub-block is within the bounds of the reference frame
            % [height, width] = size(referenceFrame);
            % refRowStart = max(1, min(refRowOffset, height - blockSize + 1));
            % refColStart = max(1, min(refColOffset, width - blockSize + 1));

            if FMEEnable
                refRowStart = 2*subBlockRowOffset -1 + mvY;
                refColStart = 2*subBlockColOffset -1 + mvX;
                referenceBlock = referenceFrame(refRowStart:2:(refRowStart + 2* blockSize - 2), refColStart:2:(refColStart + 2 * blockSize - 2));
            else
                referenceBlock = referenceFrame(refRowStart:(refRowStart + blockSize - 1), refColStart:(refColStart + blockSize - 1));
            end


            % Place the reference sub-block into the predicted block
            predictedBlock((blockY - 1) * blockSize + 1:blockY * blockSize, ...
                           (blockX - 1) * blockSize + 1:blockX * blockSize) = referenceBlock;
        end
    end
end




function [motionVector_block,total_minMAE] = compute_motionVector_block(currentBlock, currentBlockSize, referenceFrames, rowOffset, colOffset, searchRange, predictedMV, isLarge,FMEEnable, FastME)
    
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

