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
                [vector, mae, L1Norm] = findBestMatch(currentBlock, referenceFrame, row, col, blockSize, searchRange);

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
    for yOffset = -2 * searchRange:2 * searchRange
        for xOffset = -2 * searchRange:2 * searchRange
            % Determine the reference block's starting coordinates
            refRow = row + yOffset;
            refCol = col + xOffset;

            % Check if the reference block is within bounds
            if refRow > 0 && refRow + 2 * (blockSize -1) - 1 <= heightBoundry && refCol > 0 && refCol + 2 * (blockSize -1) - 1 <= widthBoundry
                % Extract the reference block
                % Note that I can only compare the size of original resolution and NOT the interpolated resolution.
                % Need to downscale the reference block by 2.
               
                %refBlock = referenceFrame(refRow:refRow + blockSize - 1, refCol:refCol + blockSize - 1);
                refBlock = referenceFrame(refRow:2:(refRow + 2* blockSize - 2), refCol:2:(refCol + 2 * blockSize - 2));
        
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