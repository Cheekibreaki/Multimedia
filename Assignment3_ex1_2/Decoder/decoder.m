
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
    params = fread(headerFile, 10, 'int32');  % [width, height, numFrames, blockSize, dct_blockSize, QP]
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
    FastME = params(10);
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
       %isIFrame = false;

        if isIFrame
            pFrameCounter = 0;  % Reset the P-frame counter
             [nonimportant1,predictionModes,quantizedResiduals,compresiduals,vbs_matrix] = entropyDecode_invquantization(isIFrame, [], encodedMDiff, encodedResidues, mvheight, mvwidth, predwidth, predheight,  reswidth, resheight,dct_blockSize, VBSEnable);
           
            
            if VBSEnable
                % compresiduals = invquantization_block(quantizedResiduals, dct_blockSize, width, height, QP,vbs_matrix);
                intraCompFrame = vbs_intraCompensation(predictionModes, compresiduals, blockSize,vbs_matrix);
            else
                % compresiduals = invquantization(quantizedResiduals, dct_blockSize, width, height, QP);
                intraCompFrame = intraCompensation(predictionModes, compresiduals, blockSize);
            end

            % Add the approximated residuals to the predicted frame to reconstruct
            reconstructedFrame = double(intraCompFrame);
            reconstructedFrame = double(max(0, min(255, reconstructedFrame)));
            interpolatedReconstructedFrame = interpolateFrame(reconstructedFrame);

            for i = 1:nRefFrames
                referenceFrames{i} = 128 * ones(height, width, 'uint8');  % Re-initialize reference frames
            end

            for i = 1:nRefFrames
                interpolatedReferenceFrames{i} = 128 * ones(height*2 - 1, width*2 -1, 'uint8');  
            end

            % if VBSEnable
            %     vbsVisualization(isIFrame,reconstructedFrame, vbs_matrix, blockSize, frameIdx, [], predictionModes);
            % else
            %     predictionInfoVisualization(reconstructedFrame, predictionModes, frameIdx, isIFrame,blockSize);
            % end
          
        else
            % [motionVectors,nonimportant1,quantizedResiduals,vbs_matrix] = entropyDecode(isIFrame, encodedMDiff, [], encodedResidues,mvheight, mvwidth, predwidth, predheight,  reswidth, resheight,dct_blockSize, VBSEnable);
             [motionVectors,nonimportant1,quantizedResiduals,compresiduals,vbs_matrix] = entropyDecode_invquantization(isIFrame, encodedMDiff, [], encodedResidues,mvheight, mvwidth, predwidth, predheight,  reswidth, resheight,dct_blockSize, VBSEnable);
            % Load the motion vectors and approximated residuals for the current frame

            if not (VBSEnable)
                motionVectors = diffDecoding(motionVectors,'mv');
            end
            % Perform motion compensation to get the predicted frame
            predictedFrame = motionCompensation(referenceFrames,interpolatedReferenceFrames, motionVectors, blockSize, width, height, FMEEnable);
            
            % if VBSEnable
            %     compresiduals = invquantization_block(quantizedResiduals, dct_blockSize, width, height, QP,vbs_matrix);
            % else
            %     compresiduals = invquantization(quantizedResiduals, dct_blockSize, width, height, QP);
            % end
            % Add the approximated residuals to the predicted frame to reconstruct
            reconstructedFrame = double(predictedFrame) + double(compresiduals);
            interpolatedReconstructedFrame = interpolateFrame(reconstructedFrame);
             
           % Save the P-frame with overlays as an image
            saveVisualizeReferenceFrames(reconstructedFrame, motionVectors, frameIdx);
            
            % if VBSEnable
            %     vbsVisualization(isIFrame,reconstructedFrame, vbs_matrix, blockSize, frameIdx, motionVectors,[]);
            % else
            %     predictionInfoVisualization(reconstructedFrame, motionVectors, frameIdx, isIFrame, blockSize);
            % end
            % 
            % Increment the P-frame counter
            pFrameCounter = min(pFrameCounter + 1, nRefFrames); 
        end
        
        % Clip the values to be in the range [0, 255] and convert to uint8
         reconstructedFrame = double(max(0, min(255, reconstructedFrame)));
        

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



function vbsVisualization(isIFrame,frame, vbs_matrix, blockSize, frameIdx, motionVectors,predictionModes)
    % Create a new figure with specified size for high resolution
    figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);  % Full-screen figure

    imshow(uint8(frame), []);
    hold on;

 % The size of the motion vector grid
    [rows, cols] = size(vbs_matrix);
    % Loop through each 2x2 block in the vbs_matrix
    for row = 1:2:rows
        for col = 1:2:cols
                % Check the 2x2 block in the vbs_matrix
                if all(vbs_matrix(row:row+1, col:col+1) == 0)
                    % Draw one large block if all 4 values are 0
                    x = ((col + 1) / 2 - 1 ) * 2 * blockSize;
                    y = ((row + 1) / 2  -1) * 2 * blockSize;
                    rectangle('Position', [x, y, 2 * blockSize, 2 * blockSize], 'EdgeColor', 'k', 'LineWidth', 1);
                   
                    if isIFrame
                         mode = predictionModes(row, col);
                         x = ((col + 1) / 2 - 1 ) * 2 * blockSize;
                         y = ((row + 1) / 2  -1) * 2 * blockSize;

                        if mode == 0
                            % Horizontal arrow
                            quiver(x , y + blockSize, 2 * blockSize, 0, 'Color', 'r', 'MaxHeadSize', 1);
                        elseif mode == 1
                            % Vertical arrow
                            quiver(x + blockSize, y, 0, 2 * blockSize, 'Color', 'b', 'MaxHeadSize', 1);
                        end
                  
                    else
                    % Draw one motion vector for the large block
                    mvRow = row;
                    mvCol = col;
                    drawMotionVector(motionVectors, mvRow, mvCol, x + blockSize, y + blockSize);
                    end
    
                else
                    x = (col - 1) * blockSize;
                    y = (row - 1) * blockSize;
                    rectangle('Position', [x, y, blockSize, blockSize], 'EdgeColor', 'k', 'LineWidth', 1);
                    rectangle('Position', [x + blockSize, y, blockSize, blockSize],  'LineWidth', 1);
                    rectangle('Position', [x, y + blockSize, blockSize, blockSize],  'LineWidth', 1);
                    rectangle('Position', [x + blockSize, y + blockSize, blockSize, blockSize],  'LineWidth', 1);
                    
                    if isIFrame
                        modeRow1 = row;
                        modeCol1 = col;
                        mode = predictionModes(modeRow1, modeCol1);
                        y = (modeRow1-1)*blockSize;
                        x = (modeCol1-1)*blockSize;
                        if mode == 0
                            % Horizontal arrow
                            quiver(x , y + blockSize/2, blockSize, 0, 'Color', 'r', 'MaxHeadSize', 1);
                        elseif mode == 1
                            % Vertical arrow
                            quiver(x + blockSize/2, y, 0, blockSize, 'Color', 'b', 'MaxHeadSize', 1);
                        end
                  
                        modeRow2 = modeRow1;
                        modeCol2 = modeCol1 +1;
                        mode = predictionModes(modeRow2, modeCol2);
                        y = (modeRow2-1)*blockSize;
                        x = (modeCol2-1)*blockSize;
                        if mode == 0
                            % Horizontal arrow
                            quiver(x , y + blockSize/2, blockSize, 0, 'Color', 'r', 'MaxHeadSize', 1);
                        elseif mode == 1
                            % Vertical arrow
                            quiver(x + blockSize/2, y, 0, blockSize, 'Color', 'b', 'MaxHeadSize', 1);
                        end

                        modeRow3 = modeRow1+1;
                        modeCol3 = modeCol1;
                        mode = predictionModes(modeRow3, modeCol3);
                        y = (modeRow3-1)*blockSize;
                        x = (modeCol3-1)*blockSize;
                        if mode == 0
                            % Horizontal arrow
                            quiver(x , y + blockSize/2, blockSize, 0, 'Color', 'r', 'MaxHeadSize', 1);
                        elseif mode == 1
                            % Vertical arrow
                            quiver(x + blockSize/2, y, 0, blockSize, 'Color', 'b', 'MaxHeadSize', 1);
                        end

                        modeRow4 = modeRow1+1;
                        modeCol4 = modeCol1+1;
                        mode = predictionModes(modeRow4, modeCol4);
                        y = (modeRow4-1)*blockSize;
                        x = (modeCol4-1)*blockSize;
                        if mode == 0
                            % Horizontal arrow
                            quiver(x , y + blockSize/2, blockSize, 0, 'Color', 'r', 'MaxHeadSize', 1);
                        elseif mode == 1
                            % Vertical arrow
                            quiver(x + blockSize/2, y, 0, blockSize, 'Color', 'b', 'MaxHeadSize', 1);
                        end

                    else
                     % Draw four motion vectors for the small blocks
                    mvRow1 = row;
                    mvCol1 = col;
                    drawMotionVector(motionVectors, mvRow1, mvCol1, x + blockSize/2, y + blockSize/2);
    
                    mvRow2 = mvRow1;
                    mvCol2 = mvCol1 + 1;
                    drawMotionVector(motionVectors, mvRow2, mvCol2, x + 3 * blockSize/2, y + blockSize/2);
    
                    mvRow3 = mvRow1 + 1;
                    mvCol3 = mvCol1;
                    drawMotionVector(motionVectors, mvRow3, mvCol3, x + blockSize/2, y + 3 * blockSize/2);
    
                    mvRow4 = mvRow1 + 1;
                    mvCol4 = mvCol1 + 1;
                    drawMotionVector(motionVectors, mvRow4, mvCol4, x + 3 * blockSize/2, y + 3 * blockSize/2);
                    end

                end
           
        end
    end

    hold off;
    % Set figure properties for high resolution
    set(gcf, 'PaperPosition', [0 0 10 10]);  % Increase the paper size
    set(gcf, 'PaperSize', [10 10]);  % Set the paper size
    set(gca, 'FontSize', 12);  % Adjust font size if needed

    % Save the visualization as a high-resolution image
    print(gcf, sprintf('../Outputs/Frame_%d_VariableBlockSize_WithPredictionInfo.png', frameIdx), '-dpng', '-r300');
    close;
end

    function predictionInfoVisualization(frame, info, frameIdx, isIFrame, blockSize)
    % Helper function to visualize prediction info
    % For I-frames, info = predictionModes (0 for horizontal, 1 for vertical)
    % For P-frames, info = motionVectors (dy, dx for each block)
     % Adjust these settings to ensure uniform arrowhead sizes and thinner lines
    arrowLineWidth = 0.5;  % Thinner lines
    
    % Create a new figure with specified size for high resolution
    figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);  % Full-screen figure

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
               
                    dy = info(by, bx, 1);
                    dx = info(by, bx, 2);
                
                x = (bx - 1) * blockSize + halfBlockSize - dx;
                y = (by - 1) * blockSize + halfBlockSize - dy;
                
                quiver(x, y, dx, dy, 'Color', 'b', ...
                      'MaxHeadSize', 1,  'LineWidth', arrowLineWidth);

            end
        end
    end

    hold off;
   % Set figure properties for high resolution
    set(gcf, 'PaperPosition', [0 0 10 10]);  % Increase the paper size
    set(gcf, 'PaperSize', [10 10]);  % Set the paper size
    set(gca, 'FontSize', 12);  % Adjust font size if needed

    % Save the visualization as a high-resolution image
    print(gcf, sprintf('../Outputs/Frame_%d_PredictionInfo.png', frameIdx), '-dpng', '-r300');
    close;
end

function drawMotionVector(motionVectors, mvRow, mvCol, endX, endY)
    % Draw a motion vector as an arrow
    dx = motionVectors(mvRow, mvCol, 1);
    dy = motionVectors(mvRow, mvCol, 2);

    startX = endX - dx;
    startY = endY - dy;

    % Scale the motion vector for better visualization
    scale = 1;  % Adjust this scale factor as needed
    quiver(startX, startY, scale * dx, scale * dy, 'Color', 'b', 'MaxHeadSize', 1, 'LineWidth', 0.5);
    
end
