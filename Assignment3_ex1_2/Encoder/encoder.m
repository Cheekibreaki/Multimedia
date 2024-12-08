function encoder(referenceFile, paddedOutputFile, numFrames, width, height, blockSize, searchRange, dct_blockSize, QP, I_Period, nRefFrames,lambda,VBSEnable,FMEEnable,FastME,RCflag,per_block_row_budget,p_bitCountPerRow,i_bitCountPerRow)
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
    lambda = lambda
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
    
    RCthreshold_bit = 100000000000000000;
    pass = 1;
    frameIdx = 1;
    over_threshold = false;
    total_per_row_qp = [];
    pass_1_motion_vectors = []
    while frameIdx <= numFrames
    
        
        
        fprintf('Processed frame %d pass %d\n', frameIdx,pass) % %d for integers
        if(RCflag > 1 && pass == 1) || RCflag ~= 2
            currentFrame = fread(fid,[width, height], 'uint8')';
        end
        %determin the frame type
        if(~over_threshold)

            if frameIdx == 1
                isIFrame = true;
            elseif frameIdx < I_Period
                isIFrame = false;
            else
                isIFrame = (mod(frameIdx - 1, I_Period) == 0);
            end
            
        else
            isIFrame = true;
            over_threshold = false;
        end
        
        %isIFrame = false;
        if isIFrame
           pFrameCounter = 0;
           if VBSEnable
               [predictedFrame, currPredictionModes, vbs_matrix,residualFrame,encodedResidues,compresiduals,total_bits_used,total_per_row_bits_used] = vbs_intraPrediction(currentFrame, blockSize,dct_blockSize,QP,lambda,RCflag,per_block_row_budget, i_bitCountPerRow,pass,total_per_row_qp);
               MDiffModes = currPredictionModes;
              
           else
               [predictedFrame, currPredictionModes] = intraPrediction(currentFrame, blockSize,dct_blockSize,QP,RCflag,per_block_row_budget, i_bitCountPerRow);
               MDiffModes = diffEncoding(currPredictionModes,'modes');
               Residuals = double(currentFrame) - double(predictedFrame);
           end
           
           
          
        else

            % Inter-frame encoding with motion estimation using multiple reference frames
            % Only use valid reference frames based on the pFrameCounter
            validRefFrames = referenceFrames(1:min(pFrameCounter + 1, nRefFrames));
            validInterpolatedRefFrames = interpolatedReferenceFrames(1:min(pFrameCounter + 1, nRefFrames));
            % Motion estimation
            if VBSEnable
                % Should use blockSize/2 here, same of motion compensation?
                [currMotionVectors, avgMAE,vbs_matrix] = vbs_motionEstimation(currentFrame, validRefFrames, validInterpolatedRefFrames, blockSize, searchRange, dct_blockSize, QP,lambda,FMEEnable, FastME, pass, pass_1_motion_vectors,RCflag);  
                MDiffMV = currMotionVectors;
                if (RCflag == 3 && pass == 1)
                    pass_1_motion_vectors = currMotionVectors;
                end
            else
                [currMotionVectors, avgMAE] = motionEstimation(currentFrame, validRefFrames,validInterpolatedRefFrames, blockSize, searchRange,FMEEnable, FastME);
                MDiffMV = diffEncoding(currMotionVectors,'mv');
            end
            % Motion compensation to get the predicted frame
            predictedFrame = motionCompensation(validRefFrames,validInterpolatedRefFrames, currMotionVectors, blockSize, width, height, FMEEnable);
            % Calculate residuals 
            Residuals = double(currentFrame) - double(predictedFrame);
            % Increment the P-frame counter
            pFrameCounter = min(pFrameCounter + 1, nRefFrames);
        end
        
             
        
        
        if isIFrame
           if VBSEnable
               quantizedResiduals = residualFrame;
                [nonimporatant1,encodedMDiff,nonimporatant1] = entropyEncode(isIFrame, [], MDiffModes, quantizedResiduals, vbs_matrix);
           else
               [quantizedResiduals,encodedResidues,compresiduals,total_bits_used,total_per_row_bits_used] = quantization_entropy(Residuals, dct_blockSize,width,height,QP,RCflag,per_block_row_budget, i_bitCountPerRow,pass,total_per_row_qp); 
               [nonimporatant1,encodedMDiff,nonimporatant1] = entropyEncode(isIFrame, [], MDiffModes, quantizedResiduals);
           end
            
            
             if(RCflag > 1 && pass == 1) 
                % if(total_bits_used > RCthreshold_bit)
                %     over_threshold = true;
                % end
                total_per_row_perc = calcRowPercent(total_per_row_bits_used,total_bits_used);
                total_per_row_qp = findQP(total_per_row_perc,QP);
                %total_per_row_qp = total_per_row_qp - 3;
                a = 1;
                
            end


            if(RCflag > 1 && pass == 2) || RCflag ~= 2
                save(sprintf('../Outputs/MDiff_frame_%d.mat', frameIdx), 'encodedMDiff');
                residualFile = sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx);
                save(residualFile, 'encodedResidues');
        
                % Clear all previous reference frames
                for i = 1:nRefFrames
                    referenceFrames{i} = 128 * ones(height, width, 'uint8');  
                end
                
                for i = 1:nRefFrames
                    interpolatedReferenceFrames{i} = 128 * ones(height*2 - 1, width*2 -1, 'uint8');  
                end
            end
      
        else
           
            if VBSEnable
                [quantizedResiduals,encodedResidues,compresiduals,total_bits_used,total_per_row_bits_used] = quantization_entropy(Residuals, dct_blockSize,width,height,QP,RCflag,per_block_row_budget, p_bitCountPerRow,pass,total_per_row_qp,vbs_matrix); 
                [encodedMDiff,nonimporatant1,nonimporatant2] = entropyEncode(isIFrame, MDiffMV, [], quantizedResiduals,vbs_matrix);
                
                
            else
                [quantizedResiduals,encodedResidues,compresiduals,total_bits_used,total_per_row_bits_used]  = quantization_entropy(Residuals, dct_blockSize,width,height,QP,RCflag,per_block_row_budget, p_bitCountPerRow,pass,total_per_row_qp); 
                [encodedMDiff,nonimporatant1,nonimporatant2] = entropyEncode(isIFrame, MDiffMV, [], quantizedResiduals);
            end
            if(RCflag > 1 && pass == 1) 
                total_per_row_perc = calcRowPercent(total_per_row_bits_used,total_bits_used);
                total_per_row_qp = findQP(total_per_row_perc,QP);
                if(total_bits_used > RCthreshold_bit)
                    over_threshold = true;
                    %total_per_row_qp = total_per_row_qp - 3;
                end
                
                
                
            end



            
            if(RCflag > 1 && pass == 2) || RCflag ~= 2

                motionVectorFile = sprintf('../Outputs/MDiff_frame_%d.mat', frameIdx);
                save(motionVectorFile, 'encodedMDiff');
                residualFile = sprintf('../Outputs/quantizedResiduals_frame_%d.mat', frameIdx);
                save(residualFile, 'encodedResidues');
            end
        end


        if(RCflag > 1 && pass == 2) || RCflag ~= 2

            reconstructedFrame = double(predictedFrame) + double(compresiduals);
            reconstructedFrame = double(max(0, min(255, reconstructedFrame)));
            interpolatedReconstructedFrame = interpolateFrame(reconstructedFrame);
    
            fwrite(yuvFile, reconstructedFrame', 'uint8');
        
            % Update the reference frames using a sliding window
            referenceFrames = [{reconstructedFrame}, referenceFrames(1:nRefFrames - 1)];
            interpolatedReferenceFrames = [{interpolatedReconstructedFrame}, interpolatedReferenceFrames(1:nRefFrames - 1)];
        end

        if(RCflag > 1 && pass == 1)
            pass = 2;
        elseif(RCflag > 1 && pass == 2)
            pass = 3;
        end

        if RCflag ~= 2
            frameIdx = frameIdx + 1;
        
        end

        if(RCflag > 1 && pass == 3)
            frameIdx = frameIdx +1;   
            pass = 1;
        end

    end
    

    fclose(fid);
    fclose(yuvFile);

end
