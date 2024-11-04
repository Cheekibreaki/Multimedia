function [motionVectors, avgMAE] = vbs_motionEstimation(currentFrame, referenceFrames, blockSize, searchRange, vbs)
    % Motion Estimation function that processes blocks in raster order.
    % Calls findBestMatch to find the best matching block in the reference frame.
    %
    % Parameters:
    %   currentFrame - The Y component of the current frame (grayscale image)
    %   referenceFrames - The Y component of the reference frames (grayscale images)
    %   blockSize - The size of the blocks (e.g., 8x8)
    %   searchRange - The maximum number of pixels to search in any direction for motion
    %   VBS_matrix - Matrix containing the VBS information (2 for single large block, 1 for four sub-blocks)
    %
    % Returns:
    %   motionVectors - An array of motion vectors for each block
    %   avgMAE - The average Mean Absolute Error for the entire frame

    % Get the dimensions of the frame
    [height, width] = size(currentFrame); % Example: 288x352
    VBS_matrix = vbs.VBS_matrix;
    % Initialize the motion vector array (for storing motion vectors for each block)
    numBlocksX = width / blockSize;
    numBlocksY = height / blockSize;
    motionVectors = zeros(numBlocksY, numBlocksX, 3);  % Stores motion vector for each block.
    groupedMatrix = zeros(numBlocksY, numBlocksX);
    totalMAE = 0;  % To keep track of the total MAE across all blocks
    
    % Process blocks in raster order (left to right, top to bottom)
    blockSkip = false(numBlocksY, numBlocksX);  % To track which blocks have already been processed
    for blockY = 1:numBlocksY
        for blockX = 1:numBlocksX
            % Skip this block if it has already been processed as part of a larger block
            if blockSkip(blockY, blockX)
                continue;
            end

            % Determine the size of the block based on VBS
            if VBS_matrix(blockY, blockX) == 2
                % Single large block: size is 2 * blockSize (i.e., twice the size in both dimensions)
                currentBlockSize = blockSize * 2;
                rowOffset = (blockY - 1) * blockSize + 1;
                colOffset = (blockX - 1) * blockSize + 1;
                
                % Extract the current block from the current frame
                currentBlock = currentFrame(rowOffset:rowOffset + currentBlockSize - 1, ...
                                            colOffset:colOffset + currentBlockSize - 1);

                % Initialize variables for the best match
                bestVector = [0, 0];
                bestRefIdx = 0;
                minMAE = inf;
                bestL1Norm = inf;

                % Check all reference frames to find the best match
                for refIdx = 1:length(referenceFrames)
                    referenceFrame = referenceFrames{refIdx};

                    % Find the best match within the current reference frame
                    [vector, mae, L1Norm] = findBestMatch(currentBlock, referenceFrame, rowOffset, colOffset, currentBlockSize, searchRange);

                    % Update if a better match is found
                    if mae < minMAE || (mae == minMAE && L1Norm < bestL1Norm)
                        minMAE = mae;
                        bestVector = vector;
                        bestL1Norm = L1Norm;
                        bestRefIdx = refIdx;
                    end
                end

                % Store the motion vector for this large block in all corresponding sub-blocks
                for dy = 0:1
                    for dx = 0:1
                        currentBlockY = blockY + dy;
                        currentBlockX = blockX + dx;

                        motionVectors(currentBlockY, currentBlockX, 1) = bestVector(1);
                        motionVectors(currentBlockY, currentBlockX, 2) = bestVector(2);
                        motionVectors(currentBlockY, currentBlockX, 3) = bestRefIdx - 1;
                        groupedMatrix(currentBlockY, currentBlockX) = 1;

                        % Mark these sub-blocks as processed
                        blockSkip(currentBlockY, currentBlockX) = true;
                    end
                end

                % Add the MAE of this block to the total MAE
                totalMAE = totalMAE + minMAE;

            elseif VBS_matrix(blockY, blockX) == 1
                % Single sub-block: size is blockSize
                currentBlockSize = blockSize;
                rowOffset = (blockY - 1) * blockSize + 1;
                colOffset = (blockX - 1) * blockSize + 1;

                % Extract the current block from the current frame
                currentBlock = currentFrame(rowOffset:rowOffset + currentBlockSize - 1, ...
                                            colOffset:colOffset + currentBlockSize - 1);

                % Initialize variables for the best match
                bestVector = [0, 0];
                bestRefIdx = 0;
                minMAE = inf;
                bestL1Norm = inf;

                % Check all reference frames to find the best match
                for refIdx = 1:length(referenceFrames)
                    referenceFrame = referenceFrames{refIdx};

                    % Find the best match within the current reference frame
                    [vector, mae, L1Norm] = findBestMatch(currentBlock, referenceFrame, rowOffset, colOffset, currentBlockSize, searchRange);

                    % Update if a better match is found
                    if mae < minMAE || (mae == minMAE && L1Norm < bestL1Norm)
                        minMAE = mae;
                        bestVector = vector;
                        bestL1Norm = L1Norm;
                        bestRefIdx = refIdx;
                    end
                end

                % Store the motion vector for this block
                motionVectors(blockY, blockX, 1) = bestVector(1);
                motionVectors(blockY, blockX, 2) = bestVector(2);
                motionVectors(blockY, blockX, 3) = bestRefIdx - 1;
                groupedMatrix(blockY, blockX) = 0;

                % Add the MAE of this block to the total MAE
                totalMAE = totalMAE + minMAE;
            end
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
