function [psnrValues, maeValues] = encoderEx3(filename, numFrames, width, height, blockSize, searchRange, n)
    % encoderEx3: This function performs motion estimation and motion 
    % compensation to encode a video sequence. It also returns PSNR and MAE values.
    %
    % Returns:
    %   psnrValues - Array of PSNR values for each frame
    %   maeValues - Array of average MAE values for each frame

    roundingFactor = 2^n;
    
    % Initialize arrays to store metrics
    psnrValues = zeros(1, numFrames);
    maeValues = zeros(1, numFrames);
    
    % Open files
    fid = fopen(filename, 'r');
    mvFile = fopen('../Outputs/motion_vectors.txt', 'w');
    yuvFile = fopen('../Outputs/referenceFrames.yuv', 'w');
    
    % First frame reference
    referenceFrame = 128 * ones(height, width, 'uint8');

    % Add array to store residual magnitudes
    residualMagnitudes = zeros(1, numFrames);

    for frameIdx = 1:numFrames
        % Read current frame
        currentFrame = fread(fid, [width, height], 'uint8')';

        if frameIdx > 1
            referenceFrame = reconstructedFrame;
        end
    
        % Motion estimation
        [motionVectors, avgMAE] = motionEstimation(currentFrame, referenceFrame, blockSize, searchRange, mvFile);
        maeValues(frameIdx) = avgMAE;
      
        % Motion compensation
        predictedFrame = motionCompensation(referenceFrame, motionVectors, blockSize);
    
        % Calculate and approximate residuals
        residuals = double(currentFrame) - double(predictedFrame);

        % Calculate residual magnitude (sum of absolute residuals)
        residualMagnitude = sum(abs(residuals(:)));
        residualMagnitudes(frameIdx) = residualMagnitude;
        
        % Print residual magnitude
        fprintf('Frame %d Residual Magnitude: %d\n', frameIdx, residualMagnitude);

        approximatedResiduals = round(residuals/roundingFactor) * roundingFactor;
    
        % Save motion vectors and residuals
        motionVectorFile = sprintf('../Outputs/motionVectors_frame_%d.mat', frameIdx);
        residualFile = sprintf('../Outputs/approximatedResiduals_frame_%d.mat', frameIdx);
        
        % Reconstruct frame
        reconstructedFrame = uint8(clip(double(predictedFrame) + approximatedResiduals));
        
        % Calculate PSNR for this frame
        mse = mean((double(currentFrame(:)) - double(reconstructedFrame(:))).^2);
        if mse > 0
            psnrValues(frameIdx) = 10 * log10(255^2/mse);
        else
            psnrValues(frameIdx) = Inf;
        end
        
        % Save reference frame
        fwrite(yuvFile, referenceFrame', 'uint8');
        save(motionVectorFile, 'motionVectors');
        save(residualFile, 'approximatedResiduals');
   
        fprintf('Processed frame %d, PSNR: %.2f dB, MAE: %.2f\n', ...
                frameIdx, psnrValues(frameIdx), maeValues(frameIdx));
    end

    fprintf('\nMagnitude of residuals table:\n');
    fprintf('Frame\t| Magnitude\n');
    fprintf('----------------------\n');
    for i = 1:numFrames
        fprintf('%d\t| %d\n', i, residualMagnitudes(i));
        % Create a figure with subplots to visualize each frame
        % --- Visualization ---

        % Calculate residuals before motion compensation (absolute value)
        residualBefore = abs(double(currentFrame) - double(referenceFrame));
        
        % Calculate residuals after motion compensation (absolute value)
        residualAfter = abs(double(currentFrame) - double(predictedFrame));
        
        % Create a figure with subplots to visualize each frame
        figure;
        % 
        %  % 1. Display the reference (previous) frame
        subplot(2, 3, 1);
        imshow(referenceFrame, []);
        title('Previous Frame');
        % 
        % % 2. Display the absolute residual before motion compensation
        subplot(2, 3, 2);
        imshow(residualBefore, []);
        title('Residual Before');
        % 
        % % 3. Display the absolute residual after motion compensation
        subplot(2, 3, 3);
        imshow(residualAfter, []);
        title('Residual After');
        % 
        % % 4. Display the source frame (current frame to encode)
        subplot(2, 3, 4);
        imshow(currentFrame, []);
        title('Source Frame');
        % 
        % % 5. Display the predicted frame after motion compensation
        subplot(2, 3, 5);
        imshow(predictedFrame, []);
        title('Predicted Frame');

    end
    
    % Close files
    fclose(fid);
    fclose(mvFile);
    fclose(yuvFile);
end

% Helper function to clip pixel values
function out = clip(in)
    out = max(0, min(255, in));
end
