function [bestVector, minMAE, bestL1Norm] = findBestMatchFractionalPixel(currentBlock, referenceFrame, row, col, blockSize, searchRange)
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

    % Convert the row/col index to interpolated frame row/col index
    row = 2*row - 1;
    col = 2*col - 1;
    % Iterate over the search range to find the best matching block in the reference frame
    for yOffset = -2 * searchRange:2 * searchRange
        for xOffset = -2 * searchRange:2 * searchRange
            % Determine the reference block's starting coordinates
            refRow = row + yOffset;
            refCol = col + xOffset;

            % Check if the reference block is within bounds
            if refRow > 0 && refRow + 2 * (blockSize -1) <= heightBoundry && refCol > 0 && refCol + 2 * (blockSize -1) <= widthBoundry
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
