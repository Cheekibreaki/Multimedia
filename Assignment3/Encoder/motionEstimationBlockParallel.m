function [motionVectors, avgMAE] = motionEstimationBlockParallel(currentFrame, originalReferenceFrames,interpolatedReferenceFrames, blockSize, searchRange, FMEEnable, FastME)
    % Motion Estimation function that processes blocks in raster order.
    % Calls findBestMatch to find the best matching block in the reference frame.
    %
    % Parameters:
    %   currentFrame - The Y component of the current frame (grayscale image)
    %   referenceFrame - The Y component of the reference frame (grayscale image)
    %   blockSize - The size of the blocks (e.g., 8x8)
    %   searchRange - The maximum number of pixels to search in any direction for motion
    %   mvFile - The file to save motion vectors
    %
    % Returns:
    %   motionVectors - An array of motion vectors for each block
    %   avgMAE - The average Mean Absolute Error for the entire frame

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

            % Store results in the 1D arrays
             motionVectors1D(idx, :) = [bestVector(1), bestVector(2), bestRefIdx - 1];
             maeValues(idx) = minMAE;
     end
    
 % Reshape motion vectors back to 2D
    motionVectors = reshape(motionVectors1D, numBlocksY, numBlocksX, 3);

    % Calculate average MAE
    avgMAE = mean(maeValues);

end

