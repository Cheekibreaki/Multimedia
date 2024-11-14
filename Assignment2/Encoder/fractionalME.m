function [motionVectors, avgMAE] = fractionalME(currentFrame, interpolatedReferenceFrames, blockSize, searchRange)
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
    
    % Initialize the motion vector array (for storing motion vectors for each block)
    numBlocksX = width / blockSize; 
    numBlocksY = height / blockSize; 
    motionVectors = zeros(numBlocksY, numBlocksX, 3);  % Stores motion vector for each block.

    totalMAE = 0;  % To keep track of the total MAE across all blocks
    
    % Process blocks in raster order (left to right, top to bottom)
    for row = 1:blockSize:height
        for col = 1:blockSize:width
            % Extract the current block from the current frame
            currentBlock = currentFrame(row:row+blockSize-1, col:col+blockSize-1);

            % Initialize variables for the best match
            bestVector = [0, 0];
            bestRefIdx = 0;
            minMAE = inf;
            bestL1Norm = inf;

            % Check all reference frames to find the best match
            for refIdx = 1:length(interpolatedReferenceFrames)
                referenceFrame = interpolatedReferenceFrames{refIdx};

                % Find the best match within the current reference frame
                [vector, mae, L1Norm] = findBestMatchFractionalPixel(currentBlock, referenceFrame, row, col, blockSize, searchRange);

                % Testing FME
                  % predictedMV = [1,1];
                  % [vector, mae, L1Norm] = findBestMatchFastFraction(currentBlock, referenceFrame, row, col, blockSize, searchRange, predictedMV);

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

            % Convert row and col to block indices
            blockY = (row - 1) / blockSize + 1;  
            blockX = (col - 1) / blockSize + 1;  

            % Store the motion vector for this block
            motionVectors(blockY, blockX, 1) = dy;
            motionVectors(blockY, blockX, 2) = dx;
            motionVectors(blockY, blockX, 3) = bestRefIdx - 1;
            % Add the MAE of this block to the total MAE
            totalMAE = totalMAE + minMAE;

           
            
        end
    end

    % Calculate the average MAE across all blocks
    avgMAE = totalMAE / (numBlocksX * numBlocksY);

end

