function decoderEx3(filename, numFrames, width, height, blockSize, searchRange)
    % decoderEx3: This function decodes the video sequence using the
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

    % Iterate through each frame to decode
    for frameIdx = 1:numFrames
        % Load the motion vectors and approximated residuals for the current frame
        motionVectorFile = sprintf('../Outputs/motionVectors_frame_%d.mat', frameIdx);
        residualFile = sprintf('../Outputs/approximatedResiduals_frame_%d.mat', frameIdx);
        load(motionVectorFile, 'motionVectors');
        load(residualFile, 'approximatedResiduals');

        % Perform motion compensation to get the predicted frame
        predictedFrame = motionCompensation(referenceFrame, motionVectors, blockSize);

        % Add the approximated residuals to the predicted frame to reconstruct
        reconstructedFrame = double(predictedFrame) + double(approximatedResiduals);

        % Clip the values to be in the range [0, 255] and convert to uint8
        reconstructedFrame = uint8(max(0, min(255, reconstructedFrame)));

        % Write the decoded frame to the output file
        fwrite(fid, reconstructedFrame', 'uint8');

        % Update the reference frame for the next iteration
        referenceFrame = reconstructedFrame;

        fprintf('Decoded frame %d\n', frameIdx);
    end

    % Close the output file
    fclose(fid);
end