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
        isIFrame = (frameIdx == 1 || mod(frameIdx - 1, I_Period) == 0);

        if isIFrame
           [predictedFrame, currPredictionModes] = intraPrediction(currentFrame, blockSize);
           % re-process the modes to be differential based on the last block
           predictionModes = diffEncoding(currPredictionModes,'modes');
         %  predictionModes = differentital(lastPredictionModes,currPredictionModes);
         %  lastPredictionModes = currPredictionModes;
        else
           
            % Motion estimation
            [currMotionVectors, avgMAE] = motionEstimation(currentFrame, referenceFrame, blockSize, searchRange);        
            % Motion compensation to get the predicted frame
            predictedFrame = motionCompensation(referenceFrame, currMotionVectors, blockSize);
            % re-process the motion vectors to be differential based on the last block
            motionVectors = diffEncoding(currMotionVectors,'mv');
        %    motionVectors = differentital(lastMotionVectors,currMotionVectors);
        %    lastMotionVectors = currMotionVectors;

        end
        
        % Calculate residuals 
        Residuals = double(currentFrame) - double(predictedFrame);
        quantizedResiduals = quantization(Residuals, dct_blockSize,width,height,QP);      
        
        % Check the quantized residual file
        txtFile = sprintf('../Outputs/testing_encoder_quantizedResiduals_frame_%d.txt',frameIdx);
        fid2 = fopen(txtFile, 'w');  % Open the file in write mode

        % Loop through the matrix and print each element
        [rows, cols] = size(quantizedResiduals);
        for i = 1:rows
            fprintf(fid2, '%d ', quantizedResiduals(i, :));  % Print elements of the row
            fprintf(fid2, '\n');  % Newline after each row
        end

        fclose(fid2);  % Close the file

      
        
        if isIFrame

            [nonimporatant1,encodedPredicitonModes,encodedResidues,predmodeBinLength,nonimporatant2] = entropyEncode(isIFrame, [], predictionModes, quantizedResiduals);
            
            save(sprintf('../Outputs/PredictionModes_frame_%d.mat', frameIdx), 'encodedPredicitonModes');
            %我们要concatinate吗？？？？？
            save(sprintf('../Outputs/predmodeBinLength_frame_%d.mat', frameIdx), 'predmodeBinLength');

            residualFile = sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx);
            save(residualFile, 'encodedResidues');

            
        else

            [encodedMotionVector,nonimporatant1,encodedResidues,nonimporatant2,motionVectorLength] = entropyEncode(isIFrame, motionVectors, [], quantizedResiduals);
            
            motionVectorFile = sprintf('../Outputs/motionVectors_frame_%d.mat', frameIdx);

            save(sprintf('../Outputs/motionVectorLength_frame_%d.mat', frameIdx), 'motionVectorLength');

            save(motionVectorFile, 'encodedMotionVector');
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



