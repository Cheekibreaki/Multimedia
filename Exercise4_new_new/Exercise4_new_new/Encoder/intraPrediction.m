function [predictedFrame, predictionModes] = intraPrediction(currentFrame, blockSize)
    [height, width] = size(currentFrame);
    predictedFrame = zeros(size(currentFrame), 'uint8');
    predictionModes = zeros(ceil(height/blockSize), ceil(width/blockSize), 'uint8');
    
    for y = 1:blockSize:height
        for x = 1:blockSize:width
            blockY = y:min(y+blockSize-1, height);
            blockX = x:min(x+blockSize-1, width);
            currentBlock = currentFrame(blockY, blockX);
            
            % Mode 0: Horizontal prediction
            if x > 1
                horizontalPred = repmat(currentFrame(blockY, x-1), 1, length(blockX));
            else
                horizontalPred = repmat(128, length(blockY), length(blockX));
            end
            
            % Mode 1: Vertical prediction
            if y > 1
                verticalPred = repmat(currentFrame(y-1, blockX), length(blockY), 1);
            else
                verticalPred = repmat(128, length(blockY), length(blockX));
            end
            
            % Mode 2: DC prediction
            if x > 1 && y > 1
                left_samples = double(currentFrame(blockY, x-1));
                top_samples = double(currentFrame(y-1, blockX));
                dcValue = round((sum(left_samples) + sum(top_samples)) / (length(left_samples) + length(top_samples)));
            elseif x > 1
                dcValue = round(mean(double(currentFrame(blockY, x-1))));
            elseif y > 1
                dcValue = round(mean(double(currentFrame(y-1, blockX))));
            else
                dcValue = 128;
            end
            dcPred = repmat(uint8(dcValue), length(blockY), length(blockX));
            
            % Calculate MAE for each mode
            maeHorizontal = mean(abs(double(currentBlock(:)) - double(horizontalPred(:))));
            maeVertical = mean(abs(double(currentBlock(:)) - double(verticalPred(:))));
            maeDC = mean(abs(double(currentBlock(:)) - double(dcPred(:))));
            
            % Choose best mode
            [~, bestMode] = min([maeHorizontal, maeVertical, maeDC]);
            predictionModes(ceil(y/blockSize), ceil(x/blockSize)) = bestMode - 1;
            
            % Set predicted block
            switch bestMode
                case 1
                    predictedFrame(blockY, blockX) = horizontalPred;
                case 2
                    predictedFrame(blockY, blockX) = verticalPred;
                case 3
                    predictedFrame(blockY, blockX) = dcPred;
            end
        end
    end
    
    % Log mode usage
    modeCount = histcounts(predictionModes, 0:3);
    fprintf('Intra Prediction Mode Usage: Horizontal: %d, Vertical: %d, DC: %d\n', modeCount);
end