function encoderEx3(filename, numFrames, width, height, blockSize, searchRange, dct_blockSize, QP, I_Period)
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


    
    % Open the padded Y only file
    fid = fopen(filename, 'r');

    % Open file for dumping motion vectors
    mvFile = fopen('../Outputs/motion_vectors.txt', 'w');
    yuvFile = fopen('../Outputs/referenceFrames.yuv', 'w');
    % For the first frame, use the hypothetical reconstructed frame as reference
    referenceFrame = 128 * ones(height, width,'uint8');  % height * width = 288 * 352
    
    for frameIdx = 1:numFrames
    
        currentFrame = fread(fid,[width, height], 'uint8')';

        if frameIdx > 1
            referenceFrame = reconstructedFrame;  % Use the reconstructed frame as the reference for the next frame
        end
        isIFrame = (frameIdx == 1 || mod(frameIdx - 1, I_Period) == 0);

        if isIFrame
           [predictedFrame, predictionModes] = intraPrediction(currentFrame, blockSize);
           save(sprintf('../Outputs/PredictionModes_frame_%d.mat', frameIdx), 'predictionModes');
        else
            % Motion estimation
            [motionVectors, avgMAE] = motionEstimation(currentFrame, referenceFrame, blockSize, searchRange, mvFile);        
            % Motion compensation to get the predicted frame
            predictedFrame = motionCompensation(referenceFrame, motionVectors, blockSize);
            motionVectorFile = sprintf('../Outputs/motionVectors_frame_%d.mat', frameIdx);
            save(motionVectorFile, 'motionVectors');
            % Calculate residuals 

        end
        
        residuals = double(currentFrame) - double(predictedFrame);
        quantizedResiduals = quantization(residuals, dct_blockSize,width,height,QP);  
        residualFile = sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx);
        
        fwrite(yuvFile, referenceFrame', 'uint8');
        
        save(residualFile, 'quantizedResiduals');
    
        % Reconstruct the frame at the encoder side to create a closed loop 
        % Use it as the reference frame for the next frame
        
        invquantizedResiduals = invquantization(quantizedResiduals, dct_blockSize,width,height,QP);
        reconstructedFrame = double(predictedFrame) + double(invquantizedResiduals);
   
        fprintf('Processed frame %d\n', frameIdx);
       
        
    end
    
    % Close the file
    fclose(fid);
    fclose(mvFile);

end



