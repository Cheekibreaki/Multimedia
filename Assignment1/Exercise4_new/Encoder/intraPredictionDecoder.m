function predictedFrame = intraPredictionDecoder(width, height, blockSize)
    predictedFrame = zeros(height, width, 'uint8');
    
    for row = 1:blockSize:height
        for col = 1:blockSize:width
            % Calculate actual block dimensions
            actualBlockHeight = min(blockSize, height-row+1);
            actualBlockWidth = min(blockSize, width-col+1);
            
            % Get left and top borders
            leftBorder = getPredictionBorder(predictedFrame, row, col-1, actualBlockHeight, 'vertical');
            topBorder = getPredictionBorder(predictedFrame, row-1, col, actualBlockWidth, 'horizontal');
            
            % Perform horizontal and vertical predictions
            hPred = horizontalPrediction(leftBorder, actualBlockWidth);
            vPred = verticalPrediction(topBorder, actualBlockHeight);
            
            % In a real implementation, we would receive the selected mode from the bitstream
            % Here, we'll just choose the mode with the lower variance as an approximation
            hVar = var(double(hPred(:)));
            vVar = var(double(vPred(:)));
            
            if hVar <= vVar
                predictedBlock = hPred;
            else
                predictedBlock = vPred;
            end
            
            % Store the predicted block
            predictedFrame(row:row+actualBlockHeight-1, col:col+actualBlockWidth-1) = predictedBlock;
        end
    end
end

function border = getPredictionBorder(frame, row, col, borderSize, direction)
    [frameHeight, frameWidth] = size(frame);
    if strcmp(direction, 'vertical')
        if col < 1
            border = 128 * ones(borderSize, 1, 'uint8');
        else
            border = frame(row:min(row+borderSize-1, frameHeight), col);
            if length(border) < borderSize
                border = [border; 128 * ones(borderSize - length(border), 1, 'uint8')];
            end
        end
    else % horizontal
        if row < 1
            border = 128 * ones(1, borderSize, 'uint8');
        else
            border = frame(row, col:min(col+borderSize-1, frameWidth));
            if length(border) < borderSize
                border = [border, 128 * ones(1, borderSize - length(border), 'uint8')];
            end
        end
    end
end

function pred = horizontalPrediction(leftBorder, width)
    pred = repmat(leftBorder, 1, width);
end

function pred = verticalPrediction(topBorder, height)
    pred = repmat(topBorder, height, 1);
end

function predictedFrame = motionCompensation(referenceFrame, motionVectors, blockSize)
    [height, width] = size(referenceFrame);
    predictedFrame = zeros(size(referenceFrame), 'uint8');
    
    for row = 1:blockSize:height
        for col = 1:blockSize:width
            blockY = ceil(row/blockSize);
            blockX = ceil(col/blockSize);
            
            dy = motionVectors(blockY, blockX, 1);
            dx = motionVectors(blockY, blockX, 2);
            
            refRowStart = max(1, min(row + dy, height - blockSize + 1));
            refColStart = max(1, min(col + dx, width - blockSize + 1));
            
            predictedFrame(row:row+blockSize-1, col:col+blockSize-1) = ...
                referenceFrame(refRowStart:refRowStart+blockSize-1, refColStart:refColStart+blockSize-1);
        end
    end
end

function reconstructedResiduals = reconstructResiduals(quantizedResiduals, QP, blockSize)
    [height, width] = size(quantizedResiduals);
    reconstructedResiduals = zeros(height, width);
    
    for row = 1:blockSize:height
        for col = 1:blockSize:width
            quantizedBlock = quantizedResiduals(row:row+blockSize-1, col:col+blockSize-1);
            dctBlock = quantizedBlock .* (2^QP);
            reconstructedBlock = idct2(dctBlock);
            reconstructedResiduals(row:row+blockSize-1, col:col+blockSize-1) = reconstructedBlock;
        end
    end
end
