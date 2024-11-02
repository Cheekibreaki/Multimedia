
function [total_bytes,bytes_list] = decoderEx4(filename)
    % decoderEx4: This function decodes the video sequence using the
    % approximated residuals and motion vectors generated during encoding.
    %
    % Parameters:
    %   filename    - File to dump decoded frames
    %   numFrames   - Number of frames to decode
    %   width       - Width of each frame
    %   height      - Height of each frame
    %   blockSize   - Size of the block for motion compensation
    %   dct_blockSize -Size of the DCT block
    %   QP          - Quantization Parameter


    % Read sequence parameters from header file
    headerFile = fopen('../Outputs/headerfile.mat', 'r');
    params = fread(headerFile, 7, 'int32');  % [width, height, numFrames, blockSize, dct_blockSize, QP]
    % Extract individual parameters
    width = params(1);
    height = params(2);
    numFrames = params(3);
    blockSize = params(4);
    dct_blockSize = params(5);
    QP = params(6);
    % Close the header file
    fclose(headerFile);
    total_bytes = 0;
    % Open file for dumping decoded frames
    fid = fopen(filename, 'w');

    % For the first frame, use the hypothetical reconstructed frame as reference
    referenceFrame = 128 * ones(height, width, 'uint8');  % height * width = 288 * 352
    mvwidth = ceil(width/blockSize);
    mvheight = ceil(height/blockSize);
    
    predwidth = ceil(width/blockSize);
    predheight = ceil(height/blockSize);

    reswidth = width;
    resheight = height;
    bytes_list = zeros(1, numFrames);
    % Iterate through each frame to decode
    for frameIdx = 1:numFrames

       quantizedresidualFile = sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx);
       load(quantizedresidualFile, 'encodedResidues');
       
       %extract the header of MDiff to get frame type
       MDiffFile = sprintf('../Outputs/MDiff_frame_%d.mat', frameIdx);
       load(MDiffFile, 'encodedMDiff');
       isIFrame = encodedMDiff(1)
       encodedMDiff = encodedMDiff(2:end); 
        quantizedResInfo = dir(quantizedresidualFile);
        MDiffInfo = dir(MDiffFile);
       total_bytes = total_bytes+quantizedResInfo.bytes
       total_bytes = total_bytes+MDiffInfo.bytes
       bytes_list(frameIdx) =  MDiffInfo.bytes + quantizedResInfo.bytes
       

        if isIFrame
            [nonimportant1,predictionModes,quantizedResiduals] = entropyDecode(isIFrame, [], encodedMDiff, encodedResidues, mvheight, mvwidth, predwidth, predheight,  reswidth, resheight);
        else
            [motionVectors,nonimportant1,quantizedResiduals] = entropyDecode(isIFrame, encodedMDiff, [], encodedResidues,mvheight, mvwidth,   predwidth, predheight,  reswidth, resheight);
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

