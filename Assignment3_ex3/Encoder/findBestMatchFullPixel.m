function [bestVector, minMAE, bestL1Norm] = findBestMatchFullPixel(currentBlock, referenceFrame, row, col, blockSize, searchRange)
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
