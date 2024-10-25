function decoderEx4(filename, numFrames, width, height, blockSize, searchRange, dct_blockSize, QP, I_Period)
    % decoderEx4: This function decodes the video sequence using the
    % approximated residuals and motion vectors generated during encoding.
    %
    % Parameters:
    %   numFrames   - Number of frames to decode
    %   width       - Width of each frame
    %   height      - Height of each frame
    %   blockSize   - Size of the block for motion compensation
    %   searchRange - Search range for motion estimation

    % Open file for dumping decoded frames
    fid = fopen(filename, 'w');

    % For the first frame, use the hypothetical reconstructed frame as reference
    referenceFrame = 128 * ones(height, width, 'uint8');  % height * width = 288 * 352
    mvwidth = ceil(width/blockSize);
    mvheight = ceil(height/blockSize);
    
    predwidth = ceil(width/blockSize);
    predheight = ceil(height/blockSize);

    % Initialize the motion vector array (for storing motion vectors for each block)
    %lastMotionVectors = zeros(mvheight, mvwidth, 2);
    %lastPredictionModes = int32(zeros(predheight, predwidth));
    
    reswidth = width;
    resheight = height;
    % Iterate through each frame to decode
    for frameIdx = 1:numFrames
        
        if frameIdx == 1
            isIFrame = true;
        elseif frameIdx < I_Period
            isIFrame = false;
        else
            isIFrame = (mod(frameIdx - 1, I_Period) == 0);
        end
        % isIFrame = false;
        
        if isIFrame
            predmodeBinLengthFile = sprintf('../Outputs/predmodeBinLength_frame_%d.mat', frameIdx);
            load(predmodeBinLengthFile, 'predmodeBinLength');
            
            predictionModesFile = sprintf('../Outputs/PredictionModes_frame_%d.mat', frameIdx);
            load(predictionModesFile, 'encodedPredicitonModes');
            quantizedresidualFile = sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx);
            load(quantizedresidualFile, 'encodedResidues');

            [nonimportant1,predictionModes,quantizedResiduals] = entropyDecode(isIFrame, [], encodedPredicitonModes, encodedResidues, mvheight, mvwidth, predwidth, predheight,  reswidth, resheight);
        else
            predictionModesFile = sprintf('../Outputs/motionVectorLength_frame_%d.mat', frameIdx);
            load(predictionModesFile, 'motionVectorLength');

            motionVectorFile = sprintf('../Outputs/motionVectors_frame_%d.mat', frameIdx);
            load(motionVectorFile, 'encodedMotionVector');
            quantizedresidualFile = sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx);
            load(quantizedresidualFile, 'encodedResidues');
            
            [motionVectors,nonimportant1,quantizedResiduals] = entropyDecode(isIFrame, encodedMotionVector, [], encodedResidues,mvheight, mvwidth,   predwidth, predheight,  reswidth, resheight);

        end


        if isIFrame
            predictionModes = diffDecoding(predictionModes,'modes');
            compresiduals = invquantization(quantizedResiduals, dct_blockSize, width, height, QP);
            intraCompFrame = intraCompensation(predictionModes, compresiduals, blockSize);


            % Add the approximated residuals to the predicted frame to reconstruct
            reconstructedFrame = double(intraCompFrame);
          
        else
            % Load the motion vectors and approximated residuals for the current frame
            motionVectors = diffDecoding(motionVectors,'mv');

            % Perform motion compensation to get the predicted frame
            predictedFrame = motionCompensation(referenceFrame, motionVectors, blockSize);
            compresiduals = invquantization(quantizedResiduals, dct_blockSize, width, height, QP);

            % Add the approximated residuals to the predicted frame to reconstruct
            reconstructedFrame = double(predictedFrame) + double(compresiduals);
        end
        
        % Clip the values to be in the range [0, 255] and convert to uint8
        % reconstructedFrame = double(max(0, min(255, reconstructedFrame)));
        

        % Write the decoded frame to the output file
        fwrite(fid, reconstructedFrame', 'uint8');

        % Update the reference frame for the next iteration
        referenceFrame = reconstructedFrame;

        fprintf('Decoded frame %d', frameIdx);
    end

    % Close the output file
    fclose(fid);
end
