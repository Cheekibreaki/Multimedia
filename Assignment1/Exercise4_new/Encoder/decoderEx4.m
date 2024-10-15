function decoderEx4(numFrames, width, height, blockSize, searchRange, QP, I_Period)
    % Initialize reconstructedFrames
    reconstructedFrames = zeros(height, width, numFrames, 'uint8');
referenceFrame = zeros(height, width, 'uint8');

for frameIdx = 1:numFrames
    fprintf('\n--- Processing frame %d ---\n', frameIdx);
    
    % Load quantized residuals
    load(sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx), 'quantizedResiduals');
    
    
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
    
    reconstructedResiduals = zeros(size(quantizedResiduals));
    for row = 1:blockSize:height
        for col = 1:blockSize:width
            blockY = row:min(row+blockSize-1, height);
            blockX = col:min(col+blockSize-1, width);
            quantizedBlock = quantizedResiduals(blockY, blockX);
            
            Q = createQMatrix(size(quantizedBlock), QP);
            dequantizedBlock = double(quantizedBlock) .* Q;
            reconstructedBlock = idct2(dequantizedBlock);
            
            reconstructedResiduals(blockY, blockX) = reconstructedBlock;
        end
    end
    
    % Scale residuals if necessary
    maxAbsResidual = max(abs(reconstructedResiduals(:)));
      if maxAbsResidual > 255
        scalingFactor = 255 / maxAbsResidual;
        reconstructedResiduals = reconstructedResiduals * scalingFactor;
      end
    
    reconstructedFrame = uint8(max(0, min(255, double(predictedFrame) + reconstructedResiduals)));
    referenceFrame = reconstructedFrame;
    
    % Print statistics
    fprintf('Frame %d statistics:\n', frameIdx);
    fprintf('Predicted frame - Min: %d, Max: %d, Mean: %.2f\n', ...
            min(predictedFrame(:)), max(predictedFrame(:)), mean(double(predictedFrame(:))));
    fprintf('Residuals - Min: %.2f, Max: %.2f, Mean: %.2f\n', ...
            min(reconstructedResiduals(:)), max(reconstructedResiduals(:)), mean(reconstructedResiduals(:)));
    fprintf('Reconstructed frame - Min: %d, Max: %d, Mean: %.2f\n', ...
            min(reconstructedFrame(:)), max(reconstructedFrame(:)), mean(double(reconstructedFrame(:))));
    
    % Display results (using the function provided earlier)
    displayDecodingResults(predictedFrame, reconstructedResiduals, reconstructedFrame, frameIdx);
end
    % Save reconstructed frames
    save('../Outputs/reconstructedFrames.mat', 'reconstructedFrames');
    fprintf('Decoder finished. Reconstructed frames saved.\n');
end

function adaptiveQP = getAdaptiveQP(block, QP)
    blockVariance = var(double(block(:)));
    if blockVariance > 1000
        adaptiveQP = max(0, QP - 2);  % Reduce QP for high-variance blocks
    elseif blockVariance < 100
        adaptiveQP = min(51, QP + 2);  % Increase QP for low-variance blocks
    else
        adaptiveQP = QP;
    end
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

function printFrameStatistics(predictedFrame, residuals, reconstructedFrame, frameIdx)
    fprintf('Frame %d statistics:\n', frameIdx);
    fprintf('Predicted frame - Min: %d, Max: %d, Mean: %.2f\n', ...
            min(predictedFrame(:)), max(predictedFrame(:)), mean(double(predictedFrame(:))));
    fprintf('Residuals - Min: %.2f, Max: %.2f, Mean: %.2f\n', ...
            min(residuals(:)), max(residuals(:)), mean(residuals(:)));
    fprintf('Reconstructed frame - Min: %d, Max: %d, Mean: %.2f\n', ...
            min(reconstructedFrame(:)), max(reconstructedFrame(:)), mean(double(reconstructedFrame(:))));
end

function displayDecodingResults(predictedFrame, residuals, reconstructedFrame, frameIdx)
    figure;
    
    % Display predicted frame
    subplot(2, 2, 1);
    imshow(predictedFrame, []);
    title('Predicted Frame');
    
    % Display residuals
    subplot(2, 2, 2);
    imshow(residuals, []); % [] allows auto-scaling for better visibility
    title('Residuals');
    colorbar; % Add a colorbar to show the scale of residuals
    
    % Display reconstructed frame
    subplot(2, 2, 3);
    imshow(reconstructedFrame, []);
    title('Reconstructed Frame');
    
    % Display histogram of reconstructed frame
    subplot(2, 2, 4);
    histogram(reconstructedFrame);
    title('Histogram of Reconstructed Frame');
    xlabel('Pixel Value');
    ylabel('Frequency');
    
    % Adjust the layout
    set(gcf, 'Position', get(0, 'Screensize')); % Maximize figure window
    
    % Pause to allow viewing (adjust time as needed)
    pause(1);
end