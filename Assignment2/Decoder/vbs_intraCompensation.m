function [reconstructedFrame] = vbs_intraCompensation(predictionModes, residuals, blockSize, vbs_matrix)

    [height, width] = size(residuals);
    reconstructedFrame = zeros(size(residuals), 'double');
    approximatedReconstructedFrame = zeros(size(residuals), 'double');
    numBlocksY = ceil(height / blockSize);
    numBlocksX = ceil(width / blockSize);
    
    for blockY = 1:2:numBlocksY
        for blockX = 1:2:numBlocksX
            % Determine whether to process as large block or split blocks
            if vbs_matrix(blockY, blockX) == 0
                % Process as large block
                [reconstructedFrame, approximatedReconstructedFrame] = processCompensationLargeBlock(...
                    reconstructedFrame, approximatedReconstructedFrame, predictionModes, residuals, ...
                    blockY, blockX, blockSize * 2);
            else
                % Process as split blocks
                for subBlockY = blockY:blockY+1
                    for subBlockX = blockX:blockX+1
                        if subBlockY <= numBlocksY && subBlockX <= numBlocksX
                            [reconstructedFrame, approximatedReconstructedFrame] = processCompensationBlock(...
                                reconstructedFrame, approximatedReconstructedFrame, predictionModes, residuals, ...
                                subBlockY, subBlockX, blockSize);
                        end
                    end
                end
            end
        end
    end
end

function [reconstructedFrame, approximatedReconstructedFrame] = processCompensationBlock(...
    reconstructedFrame, approximatedReconstructedFrame, predictionModes, residuals, ...
    blockY, blockX, blockSize)

    [height, width] = size(reconstructedFrame);

    rowOffset = (blockY - 1) * blockSize + 1;
    colOffset = (blockX - 1) * blockSize + 1;

    actualBlockHeight = min(blockSize, height - rowOffset + 1);
    actualBlockWidth = min(blockSize, width - colOffset + 1);

    % Determine the prediction mode for the current block
    predictionMode = predictionModes(blockY, blockX);

    % Initialize predicted block
    predictedBlock = zeros(actualBlockHeight, actualBlockWidth);

    if rowOffset == 1 && colOffset == 1
        % First block is predicted with mid-gray (128)
        predictedBlock(:) = 128;
    else
        switch predictionMode
            case 0  % Horizontal prediction
                if colOffset > 1
                    leftPixels = approximatedReconstructedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset - 1);
                    predictedBlock = repmat(leftPixels, 1, actualBlockWidth);
                else
                    predictedBlock(:) = 128;
                end
            case 1  % Vertical prediction
                if rowOffset > 1
                    topPixels = approximatedReconstructedFrame(rowOffset - 1, colOffset:colOffset+actualBlockWidth-1);
                    predictedBlock = repmat(topPixels, actualBlockHeight, 1);
                else
                    predictedBlock(:) = 128;
                end
            case 2  % Mid-gray prediction
                predictedBlock(:) = 128;
            otherwise
                % Handle invalid prediction modes
                predictedBlock(:) = 128;
        end
    end

    % Retrieve the residual block
    residualBlock = residuals(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1);

    % Reconstruct the block
    reconstructedBlock = double(predictedBlock) + double(residualBlock);

    % Clip pixel values to valid range [0, 255]
    reconstructedBlock = max(0, min(255, reconstructedBlock));

    % Update the reconstructed frame
    reconstructedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = reconstructedBlock;

    % Update the approximated reconstructed frame for future predictions
    approximatedReconstructedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = reconstructedBlock;
end

function [reconstructedFrame, approximatedReconstructedFrame] = processCompensationLargeBlock(...
    reconstructedFrame, approximatedReconstructedFrame, predictionModes, residuals, ...
    blockY, blockX, blockSize)

    [height, width] = size(reconstructedFrame);

    % Adjust blockY and blockX for large block positioning
    adjustedBlockY = (blockY - 1) / 2 + 1;
    adjustedBlockX = (blockX - 1) / 2 + 1;

    % Calculate base offsets for large blocks
    rowOffset = (adjustedBlockY - 1) * blockSize + 1;
    colOffset = (adjustedBlockX - 1) * blockSize + 1;

    actualBlockHeight = min(blockSize, height - rowOffset + 1);
    actualBlockWidth = min(blockSize, width - colOffset + 1);

    % Determine the prediction mode for the current block
    predictionMode = predictionModes(blockY, blockX);

    % Initialize predicted block
    predictedBlock = zeros(actualBlockHeight, actualBlockWidth);

    if rowOffset == 1 && colOffset == 1
        % First block is predicted with mid-gray (128)
        predictedBlock(:) = 128;
    else
        switch predictionMode
            case 0  % Horizontal prediction
                if colOffset > 1
                    leftPixels = approximatedReconstructedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset - 1);
                    predictedBlock = repmat(leftPixels, 1, actualBlockWidth);
                else
                    predictedBlock(:) = 128;
                end
            case 1  % Vertical prediction
                if rowOffset > 1
                    topPixels = approximatedReconstructedFrame(rowOffset - 1, colOffset:colOffset+actualBlockWidth-1);
                    predictedBlock = repmat(topPixels, actualBlockHeight, 1);
                else
                    predictedBlock(:) = 128;
                end
            case 2  % Mid-gray prediction
                predictedBlock(:) = 128;
            otherwise
                % Handle invalid prediction modes
                predictedBlock(:) = 128;
        end
    end

    % Retrieve the residual block
    residualBlock = residuals(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1);

    % Reconstruct the block
    reconstructedBlock = double(predictedBlock) + double(residualBlock);

    % Clip pixel values to valid range [0, 255]
    reconstructedBlock = max(0, min(255, reconstructedBlock));

    % Update the reconstructed frame
    reconstructedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = reconstructedBlock;

    % Update the approximated reconstructed frame for future predictions
    approximatedReconstructedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = reconstructedBlock;
end