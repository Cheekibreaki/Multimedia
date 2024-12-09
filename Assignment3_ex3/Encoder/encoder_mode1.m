function encoder_mode1(referenceFile, paddedOutputFile, numFrames, width, height, blockSize, searchRange, dct_blockSize, QP, I_Period, nRefFrames,lambda,VBSEnable,FMEEnable,FastME,mode)
    % encoder_mode1: This function performs mode 1 block level motion estimation in
    % parallel. diffencoding is disabled.

    lambda = lambda;
    % Parameters for decoder
    params.width = width;                
    params.height = height;                
    params.numFrames = numFrames;              
    params.blockSize = blockSize;                        
    params.dct_blockSize = dct_blockSize;           
    params.QP = QP;                   
    params.nRefFrames = nRefFrames;
    params.VBSEnable = VBSEnable;
    params.FMEEnable = FMEEnable;
    params.FastME = FastME;

    % Save the parameters to a MAT-file
    save('../Outputs/headerfile.mat', 'params');

    % Write parameters needed for decoder to header file
    headerFile = fopen('../Outputs/headerfile.mat', 'w');
    fwrite(headerFile, [width, height, numFrames, blockSize, dct_blockSize, QP, nRefFrames,VBSEnable, FMEEnable, FastME], 'int32');
    fclose(headerFile);

    % Open the padded Y only file
    fid = fopen(paddedOutputFile, 'r');

    % Open file for dumping motion vectors

    yuvFile = fopen(referenceFile, 'w');

    predictedFile = '../Outputs/predictedFrames.yuv';
    predictedYUVFile = fopen(predictedFile,'w');
    % % For the first frame, use the hypothetical reconstructed frame as reference
    % referenceFrame = 128 * ones(height, width,'uint8');  % height * width = 288 * 352

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

        if ~isIFrame
           

            % Inter-frame encoding with motion estimation using multiple reference frames
            % Only use valid reference frames based on the pFrameCounter
            validRefFrames = referenceFrames(1:min(pFrameCounter + 1, nRefFrames));
            validInterpolatedRefFrames = interpolatedReferenceFrames(1:min(pFrameCounter + 1, nRefFrames));
            % Motion estimation
            if VBSEnable
                [currMotionVectors, avgMAE,vbs_matrix] = vbs_motionEstimation_Mode1(mode,currentFrame, validRefFrames, validInterpolatedRefFrames, blockSize, searchRange, dct_blockSize, QP,lambda,FMEEnable, FastME);  
                MDiffMV = currMotionVectors;
            else
                [currMotionVectors, avgMAE] = motionEstimation_Mode1(currentFrame, validRefFrames,validInterpolatedRefFrames, blockSize, searchRange,FMEEnable, FastME);
                %diff encoding is disabled
                MDiffMV = currMotionVectors;
            end
            % Motion compensation to get the predicted frame
            predictedFrame = motionCompensation(validRefFrames,validInterpolatedRefFrames, currMotionVectors, blockSize, width, height, FMEEnable);
            % Calculate residuals 
            Residuals = double(currentFrame) - double(predictedFrame);
            % Increment the P-frame counter
            pFrameCounter = min(pFrameCounter + 1, nRefFrames);

        else
            % Intra-prediction is completely disabled, predictedFrame would be all grey. MDiffModes is empty
            pFrameCounter = 0;
            predictedFrame = 128 * ones(height, width,'uint8');
            Residuals = double(currentFrame) - double(predictedFrame);
            MDiffModes = [];
            % numBlocksX = width / blockSize;
            % numBlocksY = height / blockSize;
            % vbs_matrix = ones(numBlocksY, numBlocksX);
        end

        if ~ isIFrame

            if VBSEnable
                quantizedResiduals = quantization(Residuals, dct_blockSize,width,height,QP,vbs_matrix); 
                [encodedMDiff,~,encodedResiduals] = entropyEncode(mode,isIFrame, MDiffMV, [], quantizedResiduals,vbs_matrix);
            else
                quantizedResiduals = quantization(Residuals, dct_blockSize,width,height,QP); 
                [encodedMDiff,~,encodedResiduals] = entropyEncode(mode,isIFrame, MDiffMV, [], quantizedResiduals);
            end

            motionVectorFile = sprintf('../Outputs/MDiff_frame_%d.mat', frameIdx);
            save(motionVectorFile, 'encodedMDiff');
            residualFile = sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx);
            save(residualFile, 'encodedResiduals');
        else
            % if VBSEnable
            %     % For mode 1, intra is disabled, no prediction info is available
            %     quantizedResiduals = quantization(Residuals, dct_blockSize,width,height,QP,vbs_matrix); 
            %     [~,encodedMDiff,encodedResiduals] = entropyEncode(mode,isIFrame, [], MDiffModes, quantizedResiduals,vbs_matrix);
            % else
                quantizedResiduals = quantization(Residuals, dct_blockSize,width,height,QP); 
                [~,encodedMDiff,encodedResiduals] = entropyEncode(mode,isIFrame, [], MDiffModes, quantizedResiduals);
            % end
              
            % For mode 1, I frame encodedMDiff only contains the frameType header, no prediction information is transmitted
            save(sprintf('../Outputs/MDiff_frame_%d.mat', frameIdx), 'encodedMDiff');
            residualFile = sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx);
            save(residualFile, 'encodedResiduals');

            % Clear all previous reference frames
            for i = 1:nRefFrames
                referenceFrames{i} = 128 * ones(height, width, 'uint8');  
            end

            for i = 1:nRefFrames
                interpolatedReferenceFrames{i} = 128 * ones(height*2 - 1, width*2 -1, 'uint8');  
            end
                

        end
            
            % Reconstruct the frame at the encoder side to create a closed loop 
            % Use it as the reference frame for the next frame
    
            compresiduals = invquantization(quantizedResiduals, dct_blockSize,width,height,QP);
            if VBSEnable && ~isIFrame
                compresiduals = invquantization_block(quantizedResiduals, dct_blockSize, width, height, QP,vbs_matrix);
            end
            reconstructedFrame = double(predictedFrame) + double(compresiduals);
            reconstructedFrame = double(max(0, min(255, reconstructedFrame)));
            interpolatedReconstructedFrame = interpolateFrame(reconstructedFrame);
    
            fwrite(yuvFile, reconstructedFrame', 'uint8');
            fwrite(predictedYUVFile, predictedFrame','uint8');
    
            % Update the reference frames using a sliding window
            referenceFrames = [{reconstructedFrame}, referenceFrames(1:nRefFrames - 1)];
            interpolatedReferenceFrames = [{interpolatedReconstructedFrame}, interpolatedReferenceFrames(1:nRefFrames - 1)];
            fprintf('Processed frame %d\n', frameIdx);


    end

    % Close the file
    fclose(fid);
    fclose(yuvFile);
    fclose(predictedYUVFile);

end
