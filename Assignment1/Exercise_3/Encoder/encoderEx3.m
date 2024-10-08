function encoderEx3(filename, numFrames, width, height, blockSize, searchRange)
    % encoderEx3: This function performs motion estimation and motion 
    % compensation to encode a video sequence. It also visualizes the 
    % residuals before and after motion compensation for each frame.
    %
    % Parameters:
    %   filename    - The file containing the padded Y-only video frames
    %   numFrames   - Number of frames to process
    %   width       - Width of each frame
    %   height      - Height of each frame
    %   blockSize   - Size of the block for motion estimation
    %   searchRange - Search range for motion estimation

    n = 1; % You can set n = 1, 2, or 3 based on requirements
    roundingFactor = 2^n; % Rounding factor for approximated residual block (2^n)
    
    % Open the padded Y only file
    fid = fopen(filename, 'r');

    % Open file for dumping motion vectors
    mvFile = fopen('../Outputs/motion_vectors.txt', 'w');

    % For the first frame, use the hypothetical reconstructed frame as reference
    referenceFrame = 128 * ones(height, width,'uint8');  % height * width = 288 * 352
    
    for frameIdx = 1:numFrames
    
        currentFrame = fread(fid,[width, height], 'uint8')';

        if frameIdx > 1
            referenceFrame = reconstructedFrame;  % Use the reconstructed frame as the reference for the next frame
        end
    
        % Motion estimation
        [motionVectors, avgMAE] = motionEstimation(currentFrame, referenceFrame, blockSize, searchRange, mvFile);
      
        % Motion compensation to get the predicted frame
        predictedFrame = motionCompensation(referenceFrame, motionVectors, blockSize);
    
        % Calculate residuals 
        residuals = double(currentFrame) - double(predictedFrame);
    
        % Approximate the residual files
        % Will replace with quantization in exericise 4
        approximatedResiduals = round(residuals/roundingFactor) * roundingFactor;
    
        % % Save motion vectors and approximated residuals for this frame
        motionVectorFile = sprintf('../Outputs/motionVectors_frame_%d.mat', frameIdx);
        residualFile = sprintf('../Outputs/approximatedResiduals_frame_%d.mat', frameIdx);
        save(motionVectorFile, 'motionVectors');
        save(residualFile, 'approximatedResiduals');
    
        % Reconstruct the frame at the encoder side to create a closed loop 
        % Use it as the reference frame for the next frame
        reconstructedFrame = double(predictedFrame) + double(approximatedResiduals);
   
        fprintf('Processed frame %d\n', frameIdx);

        % --- Visualization ---

        % Calculate residuals before motion compensation (absolute value)
        residualBefore = abs(double(currentFrame) - double(referenceFrame));
        
        % Calculate residuals after motion compensation (absolute value)
        residualAfter = abs(double(currentFrame) - double(predictedFrame));
        
        % Create a figure with subplots to visualize each frame
        
        % figure;
        % 
        %  % 1. Display the reference (previous) frame
        % subplot(2, 3, 1);
        % imshow(referenceFrame, []);
        % title('Previous Frame');
        % 
        % % 2. Display the absolute residual before motion compensation
        % subplot(2, 3, 2);
        % imshow(residualBefore, []);
        % title('Residual Before');
        % 
        % % 3. Display the absolute residual after motion compensation
        % subplot(2, 3, 3);
        % imshow(residualAfter, []);
        % title('Residual After');
        % 
        % % 4. Display the source frame (current frame to encode)
        % subplot(2, 3, 4);
        % imshow(currentFrame, []);
        % title('Source Frame');
        % 
        % % 5. Display the predicted frame after motion compensation
        % subplot(2, 3, 5);
        % imshow(predictedFrame, []);
        % title('Predicted Frame');
    
       
        
    end
    
    % Close the file
    fclose(fid);
    fclose(mvFile);

end



