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
    refPredictionModes = zeros(ceil(height/blockSize), ceil(width/blockSize), 'uint8');
    refMotionVectors = zeros(ceil(height/blockSize), ceil(width/blockSize), 2);  % 2 channels for motion vectors (x, y)
    
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
            [predictedFrame, predictionModes] = intraPrediction(referenceFrame, blockSize);
            
            % Differential encoding of prediction modes
            diffPredictionModes = predictionModes - refPredictionModes;
            refPredictionModes = predictionModes;
            
            % Save differential prediction modes
            save(sprintf('../Outputs/diffPredictionModes_frame_%d.mat', frameIdx), 'diffPredictionModes');
        else
            fprintf('Processing P-frame\n');
            % Perform motion estimation
            [motionVectors, avgMAE] = motionEstimation(currentFrame, referenceFrame, blockSize, searchRange);
            
            % Perform motion compensation
            predictedFrame = motionCompensation(referenceFrame, motionVectors, blockSize);
            
            % Differential encoding of motion vectors
            diffMotionVectors = motionVectors - refMotionVectors;
            refMotionVectors = motionVectors;
            
            % Save differential motion vectors
            save(sprintf('../Outputs/diffMotionVectors_frame_%d.mat', frameIdx), 'diffMotionVectors');
        end    

        % Compute residuals between current and predicted frames
        residuals = double(currentFrame) - double(predictedFrame);
        
        % Quantization and DCT
        quantizedResiduals = zeros(size(residuals));
        for row = 1:dct_blockSize:height
            for col = 1:dct_blockSize:width
                block = residuals(row:row+dct_blockSize-1, col:col+dct_blockSize-1);
                dctBlock = dct2(double(block));
                
                % Quantization
                Q = createQMatrix(size(block), QP);
                quantizedBlock = round(dctBlock ./ Q);
                quantizedResiduals(row:row+dct_blockSize-1, col:col+dct_blockSize-1) = quantizedBlock;
            end
        end
        
        % Save quantized residuals
        save(sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx), 'quantizedResiduals');
        
        % Reconstruct the frame from residuals
        reconstructedResiduals = zeros(size(quantizedResiduals));
        for row = 1:blockSize:height
            for col = 1:blockSize:width
                quantizedBlock = quantizedResiduals(row:row+blockSize-1, col:col+blockSize-1);
                
                % Dequantization and inverse DCT
                Q = createQMatrix(size(quantizedBlock), QP);
                dctBlock = quantizedBlock .* Q;
                reconstructedBlock = idct2(dctBlock);
                reconstructedResiduals(row:row+blockSize-1, col:col+blockSize-1) = reconstructedBlock;
            end
        end
        
        % Clamp values to valid range
        maxAbsResidual = max(abs(reconstructedResiduals(:)));
        if maxAbsResidual > 255
            scalingFactor = 255 / maxAbsResidual;
            reconstructedResiduals = reconstructedResiduals * scalingFactor;
        end
        
        % Reconstruct the full frame
        reconstructedFrame = uint8(max(0, min(255, double(predictedFrame) + reconstructedResiduals)));
        
        % Update reference frame
        referenceFrame = reconstructedFrame;
        
        % Display the frames (original, reconstructed, residuals, and differences)
        displayEncoderFrame(currentFrame, reconstructedFrame, residuals, frameIdx);
        
        fprintf('Processed frame %d\n', frameIdx);
     end
    
    fclose(fid);  % Don't forget to close the file
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