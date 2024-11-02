function [motionVectors, avgMAE] = motionEstimation(currentFrame, referenceFrame, blockSize, searchRange)
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
    motionVectors = zeros(numBlocksY, numBlocksX, 2);  % Stores motion vector for each block.

    totalMAE = 0;  % To keep track of the total MAE across all blocks
    
    % Process blocks in raster order (left to right, top to bottom)
    for row = 1:blockSize:height
        for col = 1:blockSize:width
            % Extract the current block from the current frame
            currentBlock = currentFrame(row:row+blockSize-1, col:col+blockSize-1);

            % Find the best matching block in the reference frame using MAE
            [bestVector, minMAE] = findBestMatch(currentBlock, referenceFrame, row, col, blockSize, searchRange);

            dy = bestVector(1);
            dx = bestVector(2);

            % Convert row and col to block indices
            blockY = (row - 1) / blockSize + 1;  
            blockX = (col - 1) / blockSize + 1;  

            % Store the motion vector for this block
            motionVectors(blockY, blockX, 1) = dy;
            motionVectors(blockY, blockX, 2) = dx;

            % Add the MAE of this block to the total MAE
            totalMAE = totalMAE + minMAE;

           
            
        end
    end

    % Calculate the average MAE across all blocks
    avgMAE = totalMAE / (numBlocksX * numBlocksY);

end

function [bestVector, minMAE] = findBestMatch(currentBlock, referenceFrame, row, col, blockSize, searchRange)
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
