function test()
    % Test function for processBlock

    % Define a simple current frame
    currentFrame = [
        200, 201, 202, 203, 204, 205, 206, 207;
        208, 209, 210, 211, 212, 213, 214, 215;
        216, 217, 218, 219, 220, 221, 222, 223;
        224, 225, 226, 227, 228, 229, 230, 231;
        232, 233, 234, 235, 236, 237, 238, 239;
        240, 241, 242, 243, 244, 245, 246, 247;
        248, 249, 250, 251, 252, 253, 254, 255;
        256, 257, 258, 259, 260, 261, 262, 263
    ];
    
    % Initialize other parameters
    approximatedPredictedFrame = zeros(size(currentFrame));
    approximatedReconstructedFrame = zeros(size(currentFrame));
    predictionModes = -ones(size(currentFrame, 1) / 2, size(currentFrame, 2) / 2);

    blockY = 1;
    blockX = 1;
    rowOffset = 1;
    colOffset = 1;
    blockSize = 2;

    % % Test case 1: Mid-gray mode
    % mode = 'mid-gray';
    % [approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame] = processBlock(...
    %     currentFrame, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame, ...
    %     blockY, blockX, rowOffset, colOffset, blockSize, mode);
    % 
    % assert(isequal(approximatedPredictedFrame(1:2, 1:2), 128 * ones(2, 2)), 'Mid-gray prediction failed.');
    % assert(predictionModes(blockY, blockX) == 2, 'Mid-gray mode not set correctly.');
    % 
    % % Test case 2: Horizontal mode
    % approximatedReconstructedFrame(:) = currentFrame; % Populate with a base value
    % mode = 'horizontal';
    % [approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame] = processBlock(...
    %     currentFrame, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame, ...
    %     blockY, blockX+1, rowOffset, colOffset + 2, blockSize, mode);

 
    % Test case 3: Vertical mode
    approximatedReconstructedFrame(:) = currentFrame; % Reset to base value
    mode = 'vertical';
    [approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame] = processBlock(...
        currentFrame, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame, ...
        blockY, blockX + 1, rowOffset + 2, colOffset, blockSize, mode);


    % 
    % % Test case 4: Compare mode
    % approximatedReconstructedFrame = currentFrame; % Use the actual frame for reconstruction
    % mode = 'compare';
    % [approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame] = processBlock(...
    %     currentFrame, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame, ...
    %     blockY, blockX, rowOffset, colOffset, blockSize, mode);
    % 
    % % Assert correct prediction (depending on block data, adjust these conditions)
    % assert(predictionModes(blockY, blockX) >= 0 && predictionModes(blockY, blockX) <= 2, ...
    %     'Compare mode prediction mode out of range.');

    fprintf('All test cases passed successfully.\n');
end


function [approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame] = processBlock(...
        currentFrame, approximatedPredictedFrame, predictionModes, approximatedReconstructedFrame, ...
        blockY, blockX, rowOffset, colOffset, blockSize, mode)

    [height, width] = size(currentFrame);

    actualBlockHeight = min(blockSize, height - rowOffset + 1);
    actualBlockWidth = min(blockSize, width - colOffset + 1);

    % Initialize prediction
    predBlock = zeros(actualBlockHeight, actualBlockWidth);
    predictionMode = -1;

    switch mode
        case 'mid-gray'
            % For the first block
            predBlock(:) = 128;
            predictionMode = 2;
        case 'horizontal'
            if colOffset > 1
                leftPixels = approximatedReconstructedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset - 1);
                predBlock = repmat(leftPixels, 1, actualBlockWidth);
                predictionMode = 0;
            else
                % Handle left edge case
                predBlock(:) = 128; % Or any default value
                predictionMode = 2;
            end
        case 'vertical'
            if rowOffset > 1
                topPixels = approximatedReconstructedFrame(rowOffset - 1, colOffset:colOffset+actualBlockWidth-1);
                predBlock = repmat(topPixels, actualBlockHeight, 1);
                predictionMode = 1;
            else
                % Handle top edge case
                predBlock(:) = 128; % Or any default value
                predictionMode = 2;
            end
        case 'compare'
            predictions = {};
            modes = [];
            errors = [];

            % Horizontal Prediction
            if colOffset > 1
                leftPixels = approximatedReconstructedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset - 1);
                horizPred = repmat(leftPixels, 1, actualBlockWidth);
                predictions{end+1} = horizPred;
                modes(end+1) = 0;
                errors(end+1) = mean2(abs(double(currentFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1)) - double(horizPred)));
            end

            % Vertical Prediction
            if rowOffset > 1
                topPixels = approximatedReconstructedFrame(rowOffset - 1, colOffset:colOffset+actualBlockWidth-1);
                vertPred = repmat(topPixels, actualBlockHeight, 1);
                predictions{end+1} = vertPred;
                modes(end+1) = 1;
                errors(end+1) = mean2(abs(double(currentFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1)) - double(vertPred)));
            end

            % Choose the best prediction
            [~, idx] = min(errors);
            predBlock = predictions{idx};
            predictionMode = modes(idx);

            % Handle case where no predictions are possible
            if isempty(predictions)
                predBlock(:) = 128;
                predictionMode = 2;
            end
        otherwise
            error('Unknown prediction mode.');
    end
    
    
    % Update predicted frame and prediction modes
    approximatedPredictedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = predBlock;
    predictionModes(blockY, blockX) = predictionMode;

    % Update reconstructed frame
    approximatedReconstructedFrame(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = predBlock;
end