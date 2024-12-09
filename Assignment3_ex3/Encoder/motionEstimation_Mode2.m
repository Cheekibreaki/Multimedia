function [motionVectors, avgMAE] = motionEstimation_Mode2(currentFrame, originalReferenceFrames,interpolatedReferenceFrames, blockSize, searchRange, FMEEnable, FastME)
    % Motion Estimation function that processes blocks in parallel but keeping dependency between blocks.(Mode 2)
    % Similar to I Frame parallelism, P frame uses a similar approach, where thread 1 processes odd rows and thread 2 processes even rows.
    % Unlike I Frame, there is no vertical dependency, so thread 1 and 2 do not need to communicate with each other.

    % Get the dimensions of the frame
    [height, width] = size(currentFrame);

    % Select the appropriate reference frames
    if FMEEnable
        referenceFrames = interpolatedReferenceFrames;
    else
        referenceFrames = originalReferenceFrames;
    end

    % Initialize variables
    numBlocksX = width / blockSize;
    numBlocksY = height / blockSize;
   
    
    % Divide work among two threads
    spmd(2)
        localMotionVectors = zeros(numBlocksY, numBlocksX, 3);
        localTotalMAE = 0;
        % Determine which block rows this thread processes
        if spmdIndex == 1
            blockRowStart = 1; % Worker 1 processes odd block rows
        else
            blockRowStart = 2; % Worker 2 processes even block rows
        end
    
        % Process assigned block rows, skipping every other block row
        for blockRowIdx = blockRowStart:2:numBlocksY
            row = (blockRowIdx - 1) * blockSize + 1; % Convert block row index to pixel row index
            predictedMV = [0,0];
            for blockColIdx = 1:numBlocksX
                col = (blockColIdx - 1) * blockSize + 1; % Convert block column index to pixel column index
    
                % Extract the current block
                currentBlock = currentFrame(row:row+blockSize-1, col:col+blockSize-1);
    
                % Initialize variables for the best match
                bestRefIdx = 0;
                minMAE = inf;
                bestL1Norm = inf;
    
                % Search across all reference frames
                for refIdx = 1:length(referenceFrames)
                    referenceFrame = referenceFrames{refIdx};
    
                    % Perform motion estimation
                    if FMEEnable
                        if FastME
                            [vector, mae, L1Norm] = findBestMatchFastFraction(currentBlock, referenceFrame, row, col, blockSize, searchRange, predictedMV);
                        else
                            [vector, mae, L1Norm] = findBestMatchFractionalPixel(currentBlock, referenceFrame, row, col, blockSize, searchRange);
                        end
                    else
                        if FastME
                            [vector, mae, L1Norm] = findBestMatchFast(currentBlock, referenceFrame, row, col, blockSize, searchRange, predictedMV);
                        else
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
    
                predictedMV = bestVector;
                % Store the motion vector for this block
                localMotionVectors(blockRowIdx, blockColIdx, 1) = bestVector(1);
                localMotionVectors(blockRowIdx, blockColIdx, 2) = bestVector(2);
                localMotionVectors(blockRowIdx, blockColIdx, 3) = bestRefIdx - 1;
    
                % Add the MAE of this block to the local total MAE
                localTotalMAE = localTotalMAE + minMAE;
            end
        end
    end

    % Combine results from both threads
    motionVectors = localMotionVectors{1} + localMotionVectors{2};
    avgMAE = (localTotalMAE{1} + localTotalMAE{2}) / (numBlocksX * numBlocksY);
end