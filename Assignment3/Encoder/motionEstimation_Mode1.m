function [motionVectors, avgMAE] = motionEstimation_Mode1(currentFrame, originalReferenceFrames,interpolatedReferenceFrames, blockSize, searchRange, FMEEnable, FastME)
    % Motion Estimation function that processes all blocks in parallel.(Mode 1)

    % Get the dimensions of the frame
    [height, width] = size(currentFrame);%288 * 352

    % Check if fraction ME is enabled. If so, use the interpolated reference frames
    if FMEEnable
        referenceFrames = interpolatedReferenceFrames;
    else
        referenceFrames = originalReferenceFrames;
    end

    % Initialize the motion vector array (for storing motion vectors for each block)
    numBlocksX = width / blockSize; 
    numBlocksY = height / blockSize; 
    numBlocks = numBlocksX * numBlocksY;  % Total number of blocks

    % Preallocate motion vectors and MAE arrays
    motionVectors1D = zeros(numBlocks, 3);  % 1D array for motion vectors
    maeValues = zeros(1, numBlocks);        % 1D array for MAE values

    predictedMV = [0,0];
    
    % parallel processing all blocks
    parfor idx = 1:numBlocks 

        %translate block idx to row&col idx
        blockY = ceil(idx / numBlocksX);
        blockX = mod(idx -1,numBlocksX) + 1;
        % Calculate row and column in the frame
        row = (blockY - 1) * blockSize + 1;
        col = (blockX - 1) * blockSize + 1;

            % Extract the current block from the current frame
            currentBlock = currentFrame(row:row+blockSize-1, col:col+blockSize-1);

            % Initialize variables for the best match
            bestVector = [0, 0];
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
                    % Since diff encoding is disabled, previous vector is always [0,0]
                        [vector, mae, L1Norm] = findBestMatchFastFraction(currentBlock, referenceFrame, row, col, blockSize, searchRange, [0,0]);
                    else
                    % If fast ME is NOT enabled
                        [vector, mae, L1Norm] = findBestMatchFractionalPixel(currentBlock, referenceFrame, row, col, blockSize, searchRange);
                    end
                else
                % If fractional ME is NOT enabled
                    if FastME
                    % If fast ME is enabled
                    % Since diff encoding is disabled, previous vector is always [0,0]
                        [vector, mae, L1Norm] = findBestMatchFast(currentBlock, referenceFrame, row, col, blockSize, searchRange, [0,0]);
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

            % Store results in the 1D arrays
             motionVectors1D(idx, :) = [bestVector(1), bestVector(2), bestRefIdx - 1];
             maeValues(idx) = minMAE;
     end
    
 % Reshape motion vectors back to 2D
    motionVectors = reshape(motionVectors1D, numBlocksY, numBlocksX, 3);

        % Combine results back into motionVectors and vbs_matrix
    for idx = 1:numBlocks
       
        blockY = ceil(idx / numBlocksX);
        blockX = mod(idx -1,numBlocksX) + 1;

        motionVectors(blockY, blockX, :) = squeeze(motionVectors1D(idx, :, :, :));

    end

     % Calculate the average MAE across all blocks
    avgMAE = sum(maeValues) / numBlocks;

end

