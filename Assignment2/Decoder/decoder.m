
function [total_bytes,bytes_list] = decoder(filename)
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
    nRefFrames = params(7);
    % Close the header file
    fclose(headerFile);
    total_bytes = 0;
    % Open file for dumping decoded frames
    fid = fopen(filename, 'w');

    % Initialize a buffer to store reference frames
    referenceFrames = cell(1, nRefFrames);
    for i = 1:nRefFrames
        referenceFrames{i} = 128 * ones(height, width, 'uint8');  % Initialize reference frames
    end

    pFrameCounter = 0;  % Counter for the number of P-frames since the last Intra frame

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
       isIFrame = encodedMDiff(1);
       encodedMDiff = encodedMDiff(2:end); 
        quantizedResInfo = dir(quantizedresidualFile);
        MDiffInfo = dir(MDiffFile);
       total_bytes = total_bytes+quantizedResInfo.bytes;
       total_bytes = total_bytes+MDiffInfo.bytes;
       bytes_list(frameIdx) =  MDiffInfo.bytes + quantizedResInfo.bytes;
       

        if isIFrame
            pFrameCounter = 0;  % Reset the P-frame counter
            [nonimportant1,predictionModes,quantizedResiduals] = entropyDecode(isIFrame, [], encodedMDiff, encodedResidues, mvheight, mvwidth, predwidth, predheight,  reswidth, resheight);
            predictionModes = diffDecoding(predictionModes,'modes');
            compresiduals = invquantization(quantizedResiduals, dct_blockSize, width, height, QP);
            intraCompFrame = intraCompensation(predictionModes, compresiduals, blockSize);

            % Add the approximated residuals to the predicted frame to reconstruct
            reconstructedFrame = double(intraCompFrame);

            for i = 1:nRefFrames
                referenceFrames{i} = 128 * ones(height, width, 'uint8');  % Re-initialize reference frames
            end
          
        else
            [motionVectors,nonimportant1,quantizedResiduals] = entropyDecode(isIFrame, encodedMDiff, [], encodedResidues,mvheight, mvwidth,   predwidth, predheight,  reswidth, resheight);
            % Load the motion vectors and approximated residuals for the current frame
            % motionVectors = diffDecoding(motionVectors,'mv');

            % Perform motion compensation to get the predicted frame
            predictedFrame = motionCompensation(referenceFrames, motionVectors, blockSize);
            compresiduals = invquantization(quantizedResiduals, dct_blockSize, width, height, QP);

            % Add the approximated residuals to the predicted frame to reconstruct
            reconstructedFrame = double(predictedFrame) + double(compresiduals);

             
           % Save the P-frame with overlays as an image
            saveVisualizeReferenceFrames(reconstructedFrame, motionVectors, frameIdx);
          
            % Increment the P-frame counter
            pFrameCounter = min(pFrameCounter + 1, nRefFrames); 
        end
        
        % Clip the values to be in the range [0, 255] and convert to uint8
        % reconstructedFrame = double(max(0, min(255, reconstructedFrame)));
        

        % Write the decoded frame to the output file
        fwrite(fid, reconstructedFrame', 'uint8');

        % Update the reference frames using a sliding window
        referenceFrames = [{reconstructedFrame}, referenceFrames(1:nRefFrames - 1)];
        
        fprintf('Decoded frame %d\n', frameIdx);
    end

    % Close the output file
    fclose(fid);
end

function saveVisualizeReferenceFrames(frame, motionVectors, frameIdx)
    % Helper function to save the visualization of a frame with reference frame overlays
   % figure('Visible', 'off');  % Create an invisible figure
    figure;
    imshow(uint8(frame), []);
    hold on;

   % Define moderately contrasting colors for each reference frame index
    colors = [0.2, 0.8, 0.2;  % Medium green for index 1
              0.3, 0.3, 0.8;  % Medium blue for index 2
              0.8, 0.3, 0.3;  % Medium red for index 3
              0.7, 0.5, 0.7]; % Medium purple for index 4

    [numBlocksY, numBlocksX, ~] = size(motionVectors);
    blockSize = size(frame, 1) / numBlocksY;  % Assuming uniform block size

    % Loop through each block and draw a colored rectangle based on the reference frame index
    for by = 1:numBlocksY
        for bx = 1:numBlocksX
            refIdx = motionVectors(by, bx, 3);  % Get the reference frame index

            % Get the top-left position of the block
            x = (bx - 1) * blockSize;
            y = (by - 1) * blockSize;

            % Draw a colored rectangle
            rectangle('Position', [x, y, blockSize, blockSize], 'FaceColor', colors(refIdx + 1, :), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
    end

    hold off;
    saveas(gcf, sprintf('../Outputs/Frame_%d_Visualize_nReferenceFrames.png', frameIdx));  % Save the figure as an image
    %close; 
end