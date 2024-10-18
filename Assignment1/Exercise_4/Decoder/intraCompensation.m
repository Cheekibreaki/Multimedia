function [reconstructedFrame] = intraCompensation(predictionModes, residuals, blockSize)
    [height, width] = size(residuals);
    reconstructedFrame = zeros(size(residuals), 'uint8');
    
    for y = 1:blockSize:height
        for x = 1:blockSize:width
            actualBlockHeight = min(blockSize, height-y+1);
            actualBlockWidth = min(blockSize, width-x+1);
            
            % Determine the prediction mode for the current block
            mode = predictionModes(ceil(y/blockSize), ceil(x/blockSize));
            
            if y == 1 && x == 1
                % First block is predicted with mid-gray (128)
                predictedBlock = repmat(128, actualBlockHeight, actualBlockWidth);
            else
                % Horizontal prediction
                if mode == 0
                    if x > 1
                        predictedBlock = repmat(reconstructedFrame(y:y+actualBlockHeight-1, x-1), 1, actualBlockWidth);
                    else
                        predictedBlock = repmat(128, actualBlockHeight, actualBlockWidth);
                    end
                % Vertical prediction
                elseif mode == 1
                    if y > 1
                        predictedBlock = repmat(reconstructedFrame(y-1, x:x+actualBlockWidth-1), actualBlockHeight, 1);
                    else
                        predictedBlock = repmat(128, actualBlockHeight, actualBlockWidth);
                    end
                % Mid-gray prediction for the first block
                elseif mode == 2
                    predictedBlock = repmat(128, actualBlockHeight, actualBlockWidth);
                end
            end
            
            % Reconstruct the block by adding the residuals back to the predicted block
            reconstructedBlock = double(predictedBlock) + residuals(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1);
            
            % Store the reconstructed block in the output frame
            reconstructedFrame(y:y+actualBlockHeight-1, x:x+actualBlockWidth-1) = uint8(reconstructedBlock);
        end
    end
end