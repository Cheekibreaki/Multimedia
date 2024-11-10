
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
    params = fread(headerFile, 9, 'int32');  % [width, height, numFrames, blockSize, dct_blockSize, QP]
    % Extract individual parameters
    width = params(1);
    height = params(2);
    numFrames = params(3);
    blockSize = params(4);
    dct_blockSize = params(5);
    QP = params(6);
    nRefFrames = params(7);
    VBSEnable = params(8);
    FMEEnable = params(9);
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

    % Initialize a buffer to store interpolated reference frames.
    % Note that the interpolated size is not exactly 2*height and 2*width becasue we cannot interpolate the last row and last col.
    interpolatedReferenceFrames = cell(1, nRefFrames);
    for i = 1:nRefFrames
        interpolatedReferenceFrames{i} = 128 * ones(2*height - 1, 2*width - 1, 'uint8');
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
            % Clip the values to be in the range [0, 255] 
            reconstructedFrame = double(max(0, min(255, reconstructedFrame)));
            interpolatedReconstructedFrame = interpolateFrame(reconstructedFrame);
            
             %Save and visualize the Mode overlays
            saveVisualizePredictionInfo(reconstructedFrame, predictionModes, frameIdx, isIFrame,FMEEnable, blockSize);

            for i = 1:nRefFrames
                referenceFrames{i} = 128 * ones(height, width, 'uint8');  % Re-initialize reference frames
            end
          

             for i = 1:nRefFrames
                interpolatedReferenceFrames{i} = 128 * ones(height*2 - 1, width*2 -1, 'uint8');  
            end

        else
            [motionVectors,nonimportant1,quantizedResiduals] = entropyDecode(isIFrame, encodedMDiff, [], encodedResidues,mvheight, mvwidth,   predwidth, predheight,  reswidth, resheight);
            % Load the motion vectors and approximated residuals for the current frame
            motionVectors = diffDecoding(motionVectors,'mv');

            % Perform motion compensation to get the predicted frame
            if FMEEnable
                predictedFrame = motionCompensation(interpolatedReferenceFrames, motionVectors, blockSize, width, height, FMEEnable);
            else 
                predictedFrame = motionCompensation(referenceFrames, motionVectors, blockSize, width, height, FMEEnable);
            end
            
            compresiduals = invquantization(quantizedResiduals, dct_blockSize, width, height, QP);

            % Add the approximated residuals to the predicted frame to reconstruct
            reconstructedFrame = double(predictedFrame) + double(compresiduals);
            % Clip the values to be in the range [0, 255] 
            reconstructedFrame = double(max(0, min(255, reconstructedFrame)));
            interpolatedReconstructedFrame = interpolateFrame(reconstructedFrame);
             
           % Save the P-frame with reference frame number(color) overlays as an image
            %saveVisualizeReferenceFrames(reconstructedFrame, motionVectors, frameIdx)

            %Save and visualize the MV overlays
            saveVisualizePredictionInfo(reconstructedFrame, motionVectors, frameIdx, isIFrame, FMEEnable,blockSize);
          
            % Increment the P-frame counter
            pFrameCounter = min(pFrameCounter + 1, nRefFrames); 
        end
        
        

        % Write the decoded frame to the output file
        fwrite(fid, reconstructedFrame', 'uint8');

        % Update the reference frames using a sliding window
        referenceFrames = [{reconstructedFrame}, referenceFrames(1:nRefFrames - 1)];
        interpolatedReferenceFrames = [{interpolatedReconstructedFrame}, interpolatedReferenceFrames(1:nRefFrames - 1)];
        
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

function saveVisualizePredictionInfo(frame, info, frameIdx, isIFrame,FMEEnable, blockSize)
    % Helper function to visualize prediction info
    % For I-frames, info = predictionModes (0 for horizontal, 1 for vertical)
    % For P-frames, info = motionVectors (dy, dx for each block)
     % Adjust these settings to ensure uniform arrowhead sizes and thinner lines
    arrowLineWidth = 0.5;  % Thinner lines
    arrowHeadSize = 2;     % Uniform arrowhead size
    figure;
    imshow(uint8(frame), []);
    hold on;
    
    [numBlocksY, numBlocksX] = size(info(:,:,1));
    halfBlockSize = blockSize / 2;
    if isIFrame
        % Visualize prediction modes for I-frames
        for by = 1:numBlocksY
            for bx = 1:numBlocksX
                mode = info(by, bx);
                x = (bx - 1) * blockSize + halfBlockSize;
                y = (by - 1) * blockSize + halfBlockSize;
                
                if mode == 0
                    % Horizontal arrow
                    quiver(x - halfBlockSize, y, blockSize, 0, 'Color', 'r', 'MaxHeadSize', 1);
                elseif mode == 1
                    % Vertical arrow
                    quiver(x, y - halfBlockSize, 0, blockSize, 'Color', 'b', 'MaxHeadSize', 1);
                end
            end
        end
    else
        % Visualize motion vectors for P-frames
        for by = 1:numBlocksY
            for bx = 1:numBlocksX
                if FMEEnable
                    dy = info(by, bx, 1)/2;
                    dx = info(by, bx, 2)/2;
                else
                    dy = info(by, bx, 1);
                    dx = info(by, bx, 2);
                end
                x = (bx - 1) * blockSize + halfBlockSize;
                y = (by - 1) * blockSize + halfBlockSize;
                
                quiver(x, y, dx, dy, 'Color', 'k', ...
                       'LineWidth', arrowLineWidth, 'MaxHeadSize', arrowHeadSize);

            end
        end
    end

    hold off;
    % Save the figure as an image
    saveas(gcf, sprintf('../Outputs/PredictionInfo_Frame_%d.png', frameIdx));
    close;
end

