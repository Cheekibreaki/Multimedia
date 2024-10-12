function encoderEx4(filename, numFrames, width, height, blockSize, searchRange, QP, I_Period)
    % Calculate frame size
    frameSize = width * height;

    % Open the YUV file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end

    % Initialize reference frame
    referenceFrame = zeros(height, width, 'uint8');

    for frameIdx = 1:numFrames
        % Read the current frame
        Y = fread(fid, frameSize, 'uint8');
        if length(Y) < frameSize
            warning('End of file reached before reading all frames');
            break;
        end
        Y = reshape(Y, width, height)';
        currentFrame = Y;

        % Skip U and V components (assuming 4:2:0 format)
        fseek(fid, frameSize/2, 'cof');

        % Determine if this is an I-frame
    isIFrame = (frameIdx == 1 || mod(frameIdx - 1, I_Period) == 0);
    if isIFrame
        fprintf('Processing I-frame\n');
        predictedFrame = intraPrediction(referenceFrame, blockSize);
    else
        fprintf('Processing P-frame\n');
        % Load motion vectors and perform motion compensation
        load(sprintf('../Outputs/motionVectors_frame_%d.mat', frameIdx), 'motionVectors');
        predictedFrame = motionCompensation(referenceFrame, motionVectors, blockSize);
    end    
    residuals = double(currentFrame) - double(predictedFrame);
    
    quantizedResiduals = zeros(size(residuals));
    for row = 1:blockSize:height
        for col = 1:blockSize:width
            block = residuals(row:row+blockSize-1, col:col+blockSize-1);
            dctBlock = dct2(double(block));
            %adaptiveQP = getAdaptiveQP(block, QP);
            Q = createQMatrix(size(block), QP);
            quantizedBlock = round(dctBlock ./ Q);
            quantizedResiduals(row:row+blockSize-1, col:col+blockSize-1) = quantizedBlock;
        end
    end
    
    save(sprintf('quantizedResiduals_frame_%d.mat', frameIdx), 'quantizedResiduals');
    
    if ~isIFrame
        save(sprintf('motionVectors_frame_%d.mat', frameIdx), 'motionVectors');
    end
    
    reconstructedResiduals = zeros(size(quantizedResiduals));
    for row = 1:blockSize:height
        for col = 1:blockSize:width
            quantizedBlock = quantizedResiduals(row:row+blockSize-1, col:col+blockSize-1);
            adaptiveQP = getAdaptiveQP(quantizedBlock, QP);
            Q = createQMatrix(size(quantizedBlock), adaptiveQP);
            dctBlock = quantizedBlock .* Q;
            reconstructedBlock = idct2(dctBlock);
            reconstructedResiduals(row:row+blockSize-1, col:col+blockSize-1) = reconstructedBlock;
        end
    end

    maxAbsResidual = max(abs(reconstructedResiduals(:)));
    if maxAbsResidual > 255
        scalingFactor = 255 / maxAbsResidual;
        reconstructedResiduals = reconstructedResiduals * scalingFactor;
    end

    reconstructedFrame = uint8(max(0, min(255, double(predictedFrame) + reconstructedResiduals)));
    
    referenceFrame = reconstructedFrame;
    
    displayEncoderFrame(currentFrame, reconstructedFrame, residuals, frameIdx);
    
    fprintf('Processed frame %d\n', frameIdx);
    end
end

 function dctCoeffs = apply2DDCT(residualBlock)
        dctCoeffs = dct2(residualBlock);
 end

 function quantizedCoeffs = quantizeBlock(dctCoeffs, QP, i)
    [height, width] = size(dctCoeffs);
    Q = zeros(height, width);
    for x = 0:(i-1)
        for y = 0:(i-1)
            if (x+y < i-1)
                Q(x+1,y+1) = 2^QP;
            elseif (x+y == i-1)
                Q(x+1,y+1) = 2^(QP+1);
            else
                Q(x+1,y+1) = 2^(QP+2);
            end
        end
    end
    quantizedCoeffs = round(dctCoeffs ./ Q);
 end


 function rescaledCoeffs = rescaleBlock(quantizedCoeffs, QP, i)
    [height, width] = size(quantizedCoeffs);
    Q = zeros(height, width);
    for x = 0:(i-1)
        for y = 0:(i-1)
            if (x+y < i-1)
                Q(x+1,y+1) = 2^QP;
            elseif (x+y == i-1)
                Q(x+1,y+1) = 2^(QP+1);
            else
                Q(x+1,y+1) = 2^(QP+2);
            end
        end
    end
    rescaledCoeffs = quantizedCoeffs .* Q;
 end

 function reconstructedBlock = applyInverse2DDCT(dctCoeffs)
    reconstructedBlock = idct2(dctCoeffs);
 end

 function displayEncoderFrame(originalFrame, reconstructedFrame, residuals, frameIdx)
    figure(1);  % Use the same figure window for all frames
    
    % Display original frame
    subplot(2,2,1);
    imshow(originalFrame, []);
    title(sprintf('Original Frame %d', frameIdx));
    
    % Display reconstructed frame
    subplot(2,2,2);
    imshow(reconstructedFrame, []);
    title(sprintf('Reconstructed Frame %d', frameIdx));
    
    % Display residuals
    subplot(2,2,3);
    imshow(abs(residuals), []);  % Use abs() to visualize both positive and negative residuals
    title(sprintf('Residuals Frame %d', frameIdx));
    
    % Display difference between original and reconstructed
    subplot(2,2,4);
    difference = double(originalFrame) - double(reconstructedFrame);
    imshow(abs(difference), []);
    title(sprintf('Difference Frame %d', frameIdx));
    
    drawnow;  % Force MATLAB to render the figure immediately
    pause(0.1);  % Add a small pause to make the display visible
 end



 function filteredFrame = deblockingFilter(frame, blockSize)
    [height, width] = size(frame);
    filteredFrame = double(frame);
    
    % Vertical edges
    for x = blockSize:blockSize:width-1
        filteredFrame(:, x) = (filteredFrame(:, x-1) + 2*filteredFrame(:, x) + filteredFrame(:, x+1)) / 4;
    end
    
    % Horizontal edges
    for y = blockSize:blockSize:height-1
        filteredFrame(y, :) = (filteredFrame(y-1, :) + 2*filteredFrame(y, :) + filteredFrame(y+1, :)) / 4;
    end
    
    filteredFrame = uint8(filteredFrame);
 end

function [motionVectors, avgMAE] = motionEstimation(currentFrame, referenceFrame, blockSize, searchRange)
    [height, width] = size(currentFrame);
    motionVectors = zeros(ceil(height/blockSize), ceil(width/blockSize), 2);
    totalMAE = 0;
    blockCount = 0;

    paddedReference = padarray(referenceFrame, [searchRange searchRange], 'replicate', 'both');

    for y = 1:blockSize:height
        for x = 1:blockSize:width
            blockCount = blockCount + 1;
            currentBlock = currentFrame(y:min(y+blockSize-1,height), x:min(x+blockSize-1,width));
            [actualBlockHeight, actualBlockWidth] = size(currentBlock);
            
            bestMAE = inf;
            bestMV = [0, 0];
            
            for dy = -searchRange:searchRange
                for dx = -searchRange:searchRange
                    refY = y + searchRange + dy;
                    refX = x + searchRange + dx;
                    referenceBlock = paddedReference(refY:refY+actualBlockHeight-1, refX:refX+actualBlockWidth-1);
                    
                    MAE = sum(abs(double(currentBlock(:)) - double(referenceBlock(:)))) / (actualBlockHeight * actualBlockWidth);
                    
                    if MAE < bestMAE || (MAE == bestMAE && norm([dy, dx]) < norm(bestMV))
                        bestMAE = MAE;
                        bestMV = [dy, dx];
                    end
                end
            end
            
            motionVectors(ceil(y/blockSize), ceil(x/blockSize), :) = bestMV;
            totalMAE = totalMAE + bestMAE;
        end
    end

    avgMAE = totalMAE / blockCount;
end

function Q = createQMatrix(blockSize, QP)
    if numel(blockSize) > 1
        rows = blockSize(1);
        cols = blockSize(2);
    else
        rows = blockSize;
        cols = blockSize;
    end
    
    Q = zeros(rows, cols);
    for x = 1:rows
        for y = 1:cols
            if (x + y <= rows)
                Q(x,y) = 2^QP;
            elseif (x + y == rows + 1)
                Q(x,y) = 2^(QP+1);
            else
                Q(x,y) = 2^(QP+2);
            end
        end
    end
end

function adaptiveQP = getAdaptiveQP(block, baseQP)
    blockVariance = var(double(block(:)));
    if blockVariance > 1000
        adaptiveQP = max(0, baseQP - 2);  % Reduce QP for high-variance blocks
    elseif blockVariance < 100
        adaptiveQP = min(51, baseQP + 2);  % Increase QP for low-variance blocks
    else
        adaptiveQP = baseQP;
    end
    end