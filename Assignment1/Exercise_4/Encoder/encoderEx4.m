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
    
    % Initialize the motion vector array (for storing motion vectors for each block)

    lastMotionVectors = zeros(ceil(height/blockSize), ceil(width/blockSize),2);    
    lastQuantizedResidues = zeros(height, width);
    lastPredictionModes = int32(zeros(ceil(height/blockSize), ceil(width/blockSize)));



    for frameIdx = 1:numFrames
    
        currentFrame = fread(fid,[width, height], 'uint8')';

        if frameIdx > 1
            referenceFrame = reconstructedFrame;  % Use the reconstructed frame as the reference for the next frame
       
        end
        isIFrame = (frameIdx == 1 || mod(frameIdx - 1, I_Period) == 0);

        if isIFrame
           [predictedFrame, currPredictionModes] = intraPrediction(currentFrame, blockSize);
           predictionModes = differentital(lastPredictionModes,currPredictionModes);
           save(sprintf('../Outputs/PredictionModes_frame_%d.mat', frameIdx), 'predictionModes');
           lastPredictionModes = currPredictionModes;
        else
           
            % Motion estimation
            [currMotionVectors, avgMAE] = motionEstimation(currentFrame, referenceFrame, blockSize, searchRange, mvFile);        
            % Motion compensation to get the predicted frame
            predictedFrame = motionCompensation(referenceFrame, currMotionVectors, blockSize);
            motionVectors = differentital(lastMotionVectors,currMotionVectors);

            motionVectorFile = sprintf('../Outputs/motionVectors_frame_%d.mat', frameIdx);
            save(motionVectorFile, 'motionVectors');
            lastMotionVectors = currMotionVectors;
            % Calculate residuals 

        end
        
        residualFile = sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx);

        Residuals = double(currentFrame) - double(predictedFrame);
        currQuantizedResiduals = quantization(Residuals, dct_blockSize,width,height,QP);      
        

        quantizedResiduals = differentital(lastQuantizedResidues,currQuantizedResiduals);
        lastQuantizedResidues = currQuantizedResiduals;

        
        save(residualFile, 'quantizedResiduals');


        fwrite(yuvFile, referenceFrame', 'uint8');
        
        
    
        % Reconstruct the frame at the encoder side to create a closed loop 
        % Use it as the reference frame for the next frame
        
        invquantizedResiduals = invquantization(currQuantizedResiduals, dct_blockSize,width,height,QP);
        reconstructedFrame = double(predictedFrame) + double(invquantizedResiduals);
   
        fprintf('Processed frame %d\n', frameIdx);
       
        
    end
    
    % Close the file
    fclose(fid);
    fclose(mvFile);

end



