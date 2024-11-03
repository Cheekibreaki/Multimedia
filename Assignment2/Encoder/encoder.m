function encoder(referenceFile, paddedOutputFile, numFrames, width, height, blockSize, searchRange, dct_blockSize, QP, I_Period, nRefFrames)
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

    % Parameters for decoder
    params.width = width;                
    params.height = height;                
    params.numFrames = numFrames;              
    params.blockSize = blockSize;                        
    params.dct_blockSize = dct_blockSize;           
    params.QP = QP;                   
    params.nRefFrames = nRefFrames;
   
    % Save the parameters to a MAT-file
    save('../Outputs/headerfile.mat', 'params');

    % Write parameters needed for decoder to header file
    headerFile = fopen('../Outputs/headerfile.mat', 'w');
    fwrite(headerFile, [width, height, numFrames, blockSize, dct_blockSize, QP, nRefFrames], 'int32');
    fclose(headerFile);

    % Open the padded Y only file
    fid = fopen(paddedOutputFile, 'r');

    % Open file for dumping motion vectors

    yuvFile = fopen(referenceFile, 'w');
    
    % % For the first frame, use the hypothetical reconstructed frame as reference
    % referenceFrame = 128 * ones(height, width,'uint8');  % height * width = 288 * 352

    % Initialize a buffer to store reference frames
    referenceFrames = cell(1, nRefFrames);
    for i = 1:nRefFrames
        referenceFrames{i} = 128 * ones(height, width, 'uint8');  % Initialize reference frames
    end
    
    pFrameCounter = 0; % count number of p frames since the last intra frame. This is tracked to ensure valid number of reference frames.

    for frameIdx = 1:numFrames
    
        currentFrame = fread(fid,[width, height], 'uint8')';

        %determin the frame type
        if frameIdx == 1
            isIFrame = true;
        elseif frameIdx < I_Period
            isIFrame = false;
        else
            isIFrame = (mod(frameIdx - 1, I_Period) == 0);
        end


        if isIFrame
           pFrameCounter = 0;
           [predictedFrame, currPredictionModes] = intraPrediction(currentFrame, blockSize,dct_blockSize,QP);
           MDiffModes = diffEncoding(currPredictionModes,'modes');
           
           Residuals = double(currentFrame) - double(predictedFrame);
           
        else
            % Inter-frame encoding with motion estimation using multiple reference frames
            % Only use valid reference frames based on the pFrameCounter
            validRefFrames = referenceFrames(1:min(pFrameCounter + 1, nRefFrames));
            % Motion estimation
            [currMotionVectors, avgMAE] = motionEstimation(currentFrame, validRefFrames, blockSize, searchRange);        
            % Motion compensation to get the predicted frame
            predictedFrame = motionCompensation(validRefFrames, currMotionVectors, blockSize);
            MDiffMV = diffEncoding(currMotionVectors,'mv');
            
            % Calculate residuals 
            Residuals = double(currentFrame) - double(predictedFrame);
            % Increment the P-frame counter
            pFrameCounter = min(pFrameCounter + 1, nRefFrames);
        end
        
        
        quantizedResiduals = quantization(Residuals, dct_blockSize,width,height,QP);      
        
        
        if isIFrame

            [nonimporatant1,encodedMDiff,encodedResidues] = entropyEncode(isIFrame, [], MDiffModes, quantizedResiduals);
            
            save(sprintf('../Outputs/MDiff_frame_%d.mat', frameIdx), 'encodedMDiff');
            residualFile = sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx);
            save(residualFile, 'encodedResidues');

            % Clear all previous reference frames
            referenceFrames = cell(1, nRefFrames);
            for i = 1:nRefFrames
                referenceFrames{i} = 128 * ones(height, width, 'uint8');  
            end

      
        else

            [encodedMDiff,nonimporatant1,encodedResidues] = entropyEncode(isIFrame, MDiffMV, [], quantizedResiduals);
            
            motionVectorFile = sprintf('../Outputs/MDiff_frame_%d.mat', frameIdx);
            save(motionVectorFile, 'encodedMDiff');
            residualFile = sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx);
            save(residualFile, 'encodedResidues');
        end


        % Reconstruct the frame at the encoder side to create a closed loop 
        % Use it as the reference frame for the next frame
        
        compresiduals = invquantization(quantizedResiduals, dct_blockSize,width,height,QP);
        reconstructedFrame = double(predictedFrame) + double(compresiduals);

        fwrite(yuvFile, reconstructedFrame', 'uint8');

        % Update the reference frames using a sliding window
        referenceFrames = [{reconstructedFrame}, referenceFrames(1:nRefFrames - 1)];

        fprintf('Processed frame %d\n', frameIdx);
       
        
    end
    
    % Close the file
    fclose(fid);
    fclose(yuvFile);

end
