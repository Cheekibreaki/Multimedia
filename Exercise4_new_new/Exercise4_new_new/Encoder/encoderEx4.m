function encoderEx4(filename, numFrames, width, height, blockSize, dct_blockSize, searchRange, QP, I_Period)
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
        %predictedFrame = intraPrediction(referenceFrame, blockSize);
        [predictedFrame, predictionModes] = intraPrediction(currentFrame, blockSize);
        % Save prediction modes
        save(sprintf('../Outputs/predictionModes_frame_%d.mat', frameIdx), 'predictionModes');
    else
        fprintf('Processing P-frame\n');
        % Load motion vectors and perform motion compensation
        load(sprintf('../Outputs/motionVectors_frame_%d.mat', frameIdx), 'motionVectors');
        predictedFrame = motionCompensation(referenceFrame, motionVectors, blockSize);
    end    
    residuals = double(currentFrame) - double(predictedFrame);
    
    quantizedResiduals = zeros(size(residuals));
    for row = 1:dct_blockSize:height
        for col = 1:dct_blockSize:width
            block = residuals(row:row+dct_blockSize-1, col:col+dct_blockSize-1);
            dctBlock = dct2(double(block));
            %debug
            %fprintf('DCT Coefficients - Min: %.2f, Max: %.2f, Mean: %.2f\n', min(dctBlock(:)), max(dctBlock(:)), mean(dctBlock(:)));
            Q = createQMatrix(size(block), QP);
            quantizedBlock = round(dctBlock ./ Q);
            %quantizedResiduals(row:row+blockSize-1, col:col+blockSize-1) = quantizedBlock;
            quantizedResiduals(row:min(row+dct_blockSize-1,height), col:min(col+dct_blockSize-1,width)) = quantizedBlock;
        end
    end
    
    % Save quantized residuals
    save(sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx), 'quantizedResiduals');
    
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
    fprintf('Quantized Coefficients - Min: %d, Max: %d, Mean: %.2f\n', min(quantizedCoeffs(:)), max(quantizedCoeffs(:)), mean(quantizedCoeffs(:)));
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
    fprintf('Dequantized Coefficients - Min: %.2f, Max: %.2f, Mean: %.2f\n', min(dequantizedCoeffs(:)), max(dequantizedCoeffs(:)), mean(dequantizedCoeffs(:)));
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









