function encoderEx4(referenceFile, paddedOutputFile, numFrames, width, height, blockSize, searchRange, dct_blockSize, QP, I_Period)
    % encoderEx3: This function performs motion estimation and motion 
    % compensation to encode a video sequence. It also visualizes the 
    % residuals before and after motion compensation for each frame.
    %
    % Parameters:
    %   paddedOutputFile    - The file containing the padded Y-only video frames
    %   numFrames   - Number of frames to process
    %   width       - Width of each frame
    %   height      - Height of each frame
    %   blockSize   - Size of the block for motion estimation
    %   searchRange - Search range for motion estimation


    
    % Open the padded Y only file
    fid = fopen(paddedOutputFile, 'r');

    % Open file for dumping motion vectors

    yuvFile = fopen(referenceFile, 'w');
    
    % For the first frame, use the hypothetical reconstructed frame as reference
    referenceFrame = 128 * ones(height, width,'uint8');  % height * width = 288 * 352
    
    % Initialize the motion vector array (for storing motion vectors for each block)

    %lastMotionVectors = zeros(ceil(height/blockSize), ceil(width/blockSize),2);    
    %lastPredictionModes = int32(zeros(ceil(height/blockSize), ceil(width/blockSize)));



    for frameIdx = 1:numFrames
    
        currentFrame = fread(fid,[width, height], 'uint8')';

        if frameIdx > 1
            referenceFrame = reconstructedFrame;  % Use the reconstructed frame as the reference for the next frame
       
        end

        %determin the frame type
        if frameIdx == 1
            isIFrame = true;
        elseif frameIdx < I_Period
            isIFrame = false;
        else
            isIFrame = (mod(frameIdx - 1, I_Period) == 0);
        end


        if isIFrame
           [predictedFrame, currPredictionModes] = intraPrediction(currentFrame, blockSize,dct_blockSize,QP);
           MDiffModes = diffEncoding(currPredictionModes,'modes');
           
           Residuals = double(currentFrame) - double(predictedFrame);
           
        else
           
            % Motion estimation
            [currMotionVectors, avgMAE] = motionEstimation(currentFrame, referenceFrame, blockSize, searchRange);        
            % Motion compensation to get the predicted frame
            predictedFrame = motionCompensation(referenceFrame, currMotionVectors, blockSize);
            MDiffMV = diffEncoding(currMotionVectors,'mv');
            
            % Calculate residuals 
            Residuals = double(currentFrame) - double(predictedFrame);
        end
        
        
        
        quantizedResiduals = quantization(Residuals, dct_blockSize,width,height,QP);      
        
        
        if isIFrame

            [nonimporatant1,encodedMDiff,encodedResidues] = entropyEncode(isIFrame, [], MDiffModes, quantizedResiduals);
            
            save(sprintf('../Outputs/MDiff_frame_%d.mat', frameIdx), 'encodedMDiff');
            residualFile = sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx);
            save(residualFile, 'encodedResidues');

      
        else

            [encodedMDiff,nonimporatant1,encodedResidues] = entropyEncode(isIFrame, MDiffMV, [], quantizedResiduals);
            
            motionVectorFile = sprintf('../Outputs/MDiff_frame_%d.mat', frameIdx);
            save(motionVectorFile, 'encodedMDiff');
            residualFile = sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx);
            save(residualFile, 'encodedResidues');
        end

        
        fwrite(yuvFile, referenceFrame', 'uint8');
        
        % Reconstruct the frame at the encoder side to create a closed loop 
        % Use it as the reference frame for the next frame
        
        compresiduals = invquantization(quantizedResiduals, dct_blockSize,width,height,QP);
        reconstructedFrame = double(predictedFrame) + double(compresiduals);
        fprintf('Processed frame %d\n', frameIdx);
       
        
    end
    
    % Close the file
    fclose(fid);

end



