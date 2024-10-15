function predictedFrame = intraCompensation(predictionModes, blockSize, height, width)
    [modeHeight, modeWidth] = size(predictionModes);
    predictedFrame = zeros(height, width, 'uint8');

   % SOMETHING IS WRONG HERE, it shouldn't be zeros, should be last ref
   % frame
    
    for y = 1:modeHeight
        for x = 1:modeWidth
            mode = predictionModes(y, x);
            block = predictBlock(predictedFrame, x, y, blockSize, mode);
            yStart = (y-1)*blockSize + 1;
            xStart = (x-1)*blockSize + 1;
            yEnd = min(y*blockSize, height);
            xEnd = min(x*blockSize, width);
            predictedFrame(yStart:yEnd, xStart:xEnd) = block(1:(yEnd-yStart+1), 1:(xEnd-xStart+1));
        end
    end
end

function block = predictBlock(frame, x, y, blockSize, mode)
    [height, width] = size(frame);
    blockY = (y-1)*blockSize+1 : min(y*blockSize, height);
    blockX = (x-1)*blockSize+1 : min(x*blockSize, width);
    
    switch mode
        case 0 % Horizontal
            if (x-1)*blockSize > 0
                block = repmat(frame(blockY, (x-1)*blockSize), 1, length(blockX));
            else
                block = repmat(128, length(blockY), length(blockX));
            end
        case 1 % Vertical
            if (y-1)*blockSize > 0
                block = repmat(frame((y-1)*blockSize, blockX), length(blockY), 1);
            else
                block = repmat(128, length(blockY), length(blockX));
            end
        case 2 % DC
            if (x-1)*blockSize > 0 && (y-1)*blockSize > 0
                left_samples = double(frame(blockY, (x-1)*blockSize));
                top_samples = double(frame((y-1)*blockSize, blockX));
                dc_val = round((sum(left_samples) + sum(top_samples)) / (length(left_samples) + length(top_samples)));
            elseif (x-1)*blockSize > 0
                dc_val = round(mean(double(frame(blockY, (x-1)*blockSize))));
            elseif (y-1)*blockSize > 0
                dc_val = round(mean(double(frame((y-1)*blockSize, blockX))));
            else
                dc_val = 128;
            end
            block = repmat(uint8(dc_val), length(blockY), length(blockX));
        otherwise
            error('Invalid prediction mode: %d', mode);
    end
end