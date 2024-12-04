function [bestVector, minMAE, bestL1Norm] = findBestMatchFast(currentBlock, referenceFrame, row, col, blockSize, searchRange, predictedMV)
    % Nearest Neighbors Search for Motion Estimation.
    % 
    % Parameters:
    %   currentBlock - The block from the current frame to match
    %   referenceFrame - The reference frame 
    %   row, col - The starting position of the block in the current frame
    %   blockSize - The size of the block
    %   searchRange - The maximum pixel offset to search in the reference frame
    %   predictedMV - The predicted motion vector (dy, dx) from the previous block
    %
    % Returns:
    %   bestVector - The best motion vector (dy, dx) for the current block
    %   minMAE - The minimum MAE for the best matching block
    %   bestL1Norm - The L1 norm of the best motion vector

    % Initialization
    minMAE = inf; 
    bestVector = [0, 0]; 
    bestL1Norm = inf;

    % Reference Frame boundaries
    [heightBoundry, widthBoundry] = size(referenceFrame);

    % Step 1: Search the (0, 0) location
    refRow = row;
    refCol = col;

    % Check if the (0, 0) location is within bounds
    if refRow > 0 && refRow + blockSize - 1 <= heightBoundry && refCol > 0 && refCol + blockSize - 1 <= widthBoundry
        % Extract the reference block at (0, 0)
        refBlock = referenceFrame(refRow:refRow + blockSize - 1, refCol:refCol + blockSize - 1);

        % Calculate the Mean Absolute Error (MAE)
        currentBlock = double(currentBlock);
        refBlock = double(refBlock);
        MAE = mean(abs(currentBlock(:) - refBlock(:)));

        % Compute the L1 norm of the motion vector (0, 0)
        L1Norm = abs(0) + abs(0);

        % Update the best match
        minMAE = MAE;
        bestVector = [0, 0];
        bestL1Norm = L1Norm;
    end

    % Step 2: Search starting from the predicted motion vector
    searchOrigin = predictedMV;
    keepSearching = true;

    % Iterate to refine the search around the best match
    while keepSearching
        % Search the origion and four neighboring positions around the
        % current origin witha '+' shape
        candidateVectors = [
            searchOrigin;
            searchOrigin + [0, 1];   % Right
            searchOrigin + [0, -1];  % Left
            searchOrigin + [1, 0];   % Down
            searchOrigin + [-1, 0]   % Up
        ];

        currentBestVector = searchOrigin;
        currentMinMAE = inf;
        currentBestL1Norm = inf;

        % Evaluate each candidate motion vector
        for k = 1:size(candidateVectors, 1)
            yOffset = candidateVectors(k, 1);
            xOffset = candidateVectors(k, 2);

            % Added condition to ensure the motion vector stays within the search range 
            if abs(yOffset) > searchRange || abs(xOffset) > searchRange
                continue; % Skip this candidate if it exceeds the allowed range
            end

            % Determine the reference block's starting coordinates
            refRow = row + yOffset;
            refCol = col + xOffset;

            % Check if the reference block is within bounds
            if refRow > 0 && refRow + blockSize - 1 <= heightBoundry && refCol > 0 && refCol + blockSize - 1 <= widthBoundry
                % Extract the reference block
                refBlock = referenceFrame(refRow:refRow + blockSize - 1, refCol:refCol + blockSize - 1);

                % Calculate the Mean Absolute Error (MAE)
                refBlock = double(refBlock);
                MAE = mean(abs(currentBlock(:) - refBlock(:)));

                % Compute the L1 norm of the current motion vector
                currentL1Norm = abs(yOffset) + abs(xOffset);

                % Update if a better match is found
                if MAE < currentMinMAE || (MAE == currentMinMAE && currentL1Norm < currentBestL1Norm)
                    currentMinMAE = MAE;
                    currentBestVector = [yOffset, xOffset];
                    currentBestL1Norm = currentL1Norm;
                end
            end
        end

        % Check if the best match is at the current origin
        if currentBestVector(1) == searchOrigin(1) && currentBestVector(2) == searchOrigin(2)
            keepSearching = false;  % Stop searching if the origin did not change
        else
            % Set the new origin to the position of the best match
            searchOrigin = currentBestVector;
        end
    end

    % Compare the best motion vector found during the nearest neighbors search with MV at (0, 0), Update the overall best match if necessary
    if currentMinMAE < minMAE || (currentMinMAE == minMAE && currentBestL1Norm < bestL1Norm)
        minMAE = currentMinMAE;
        bestVector = currentBestVector;
        bestL1Norm = currentBestL1Norm;
    end
end

