function [decodedMotionVector3d,decodedPredicitonModes2d,decodedResidues2d,reconstructedResiduals2d,reconstructed_vbs_matrix] = entropyDecode_invquantization(frame_type, encodedMotionVector1d, encodedPredicitonModes1d, encodedResidues1d,mvwidth, mvheight,  predwidth, predheight,  reswidth, resheight,dct_blockSize, VBSEnable)
    
    %switch the order of width & height
    temp = mvwidth;
    mvwidth = mvheight;
    mvheight = temp;

    temp = predwidth;
    predwidth = predheight;
    predheight = temp;
    
    temp = reswidth;
    reswidth = resheight;
    resheight = temp;

    % Input:
    % frame_type: 1 for I-frame, 0 for P-frame 
    % pred_diff: array containing differential prediction information (modes for intra or motion vectors for inter)
    % dct_coeffs: block of quantized DCT coefficients
    % block_size: size of the block (e.g., 4 for 4x4 blocks)
    reconstructed_vbs_matrix = [];  % VBS not used
    decodedMotionVector3d = [];
    decodedPredicitonModes2d = [];
    decodedResidues2d = [];
    if frame_type == 0  % Assuming data handling for P-frames
        if VBSEnable
            % Decoding motion vectors using the structure from the encoding algorithm
            motionVectorRevEGC = exp_golomb_decode(encodedMotionVector1d);
            idx = 1;
            [rows, cols] = deal(mvheight, mvwidth);
            decodedMotionVector3d = zeros(rows, cols, 3);
            previous_motion_vector_block = zeros(1,1,3);
            reconstructed_vbs_matrix = zeros(rows, cols);
            % Loop through the decoded values in a 2x2 block size
            for row = 1:2:rows
                for col = 1:2:cols
                    % Extract prefix (0 or 1) from the encoded data
                    prefix = motionVectorRevEGC(idx);
                    idx = idx + 1;
                    
                    if prefix == 1
                        % Decode the entire 2x2x3 block
                        motionVector1d = motionVectorRevEGC(idx:idx+11);
                        idx = idx + 12;
                        motionVector2d = invzigzag(motionVector1d, 2,2 * 3);
                        motion_block = reshape_2d_to_3d(motionVector2d, 2, 2, 3);
                        [motion_block,previous_motion_vector_block] = diffDecoding_block(motion_block,'mv',previous_motion_vector_block);
                        decodedMotionVector3d(row:row+1, col:col+1, :) = motion_block;
                        reconstructed_vbs_matrix(row:row+1, col:col+1) = 1;
                    elseif prefix == 0
                        % Decode the top-left 1x1x3 element
                        motionVector1d = motionVectorRevEGC(idx:idx+2);
                        idx = idx + 3;
                        motionVector2d = invzigzag(motionVector1d, 1, 3);
                        top_left_block = reshape_2d_to_3d(motionVector2d, 1, 1, 3);
                        % Assign the top-left block to all elements in the 2x2 block
                        
                        motion_block = zeros(2,2,3); 
                        motion_block(1,1,:) = top_left_block;
    
                        [motion_block,previous_motion_vector_block] = diffDecoding_block(motion_block,'mv',previous_motion_vector_block);
                        decodedMotionVector3d(row:row+1, col:col+1, :) = motion_block;
                        reconstructed_vbs_matrix(row:row+1, col:col+1) = 0;
                    else
                        error('ErrorID:1', 'Invalid prefix encountered during decoding!');
                    end
                    
                end
            end

        else
            motionVectorRevEGC = exp_golomb_decode(encodedMotionVector1d);
            motionVector2d = invzigzag(motionVectorRevEGC, mvheight, length(motionVectorRevEGC)/mvheight);
            decodedMotionVector3d = reshape_2d_to_3d(motionVector2d, mvwidth, mvheight, 3);
        end

    elseif frame_type == 1
        if VBSEnable
            temp = predwidth;
            predwidth = predheight;
            predheight = temp;
           

            % Decoding prediction modes using the structure from the encoding algorithm
            predictedModeRevEGC = exp_golomb_decode(encodedPredicitonModes1d);
            idx = 1;
            [rows, cols] = deal(predheight, predwidth);
            decodedPredicitonModes2d = zeros(rows, cols);
            reconstructed_vbs_matrix = zeros(rows, cols);
            previous_prediction_mode_block = 0;

            % Loop through the decoded values in a 2x2 block size
            for row = 1:2:rows
                for col = 1:2:cols
                    % Extract prefix (0 or 1) from the encoded data
                    if idx > length(predictedModeRevEGC)
                        break;  % Reached the end of the encoded data
                    end
                    prefix = predictedModeRevEGC(idx);
                    idx = idx + 1;

                    if prefix == 1  % Split block
                        % Read the next 4 elements (2x2 block)
                        num_elements = 4;  % Since zigzag of 2x2 block is 4 elements
                        if idx + num_elements - 1 > length(predictedModeRevEGC)
                            error('Not enough data to decode split block');
                        end
                        predictionMode1d = predictedModeRevEGC(idx:idx + num_elements -1);
                        idx = idx + num_elements;

                        % Inverse zigzag to get 2x2 block
                        predictionMode_block = invzigzag(predictionMode1d, 2, 2);

                        % Apply differential decoding
                        [predictionMode_block, previous_prediction_mode_block] = ...
                            diffDecoding_block(predictionMode_block, 'modes', previous_prediction_mode_block);

                        % Place the decoded modes into the correct positions
                        decodedPredicitonModes2d(row:row+1, col:col+1) = predictionMode_block;

                        % Set vbs_matrix to 1(split block)
                        reconstructed_vbs_matrix(row:row+1, col:col+1) = 1;

                    elseif prefix == 0  % Large block
                        % Read the next element (top-left mode)
                        if idx > length(predictedModeRevEGC)
                            error('Not enough data to decode large block');
                        end
                        predictionMode1d = predictedModeRevEGC(idx);
                        idx = idx + 1;

                        % Inverse zigzag to get the mode (only one element)
                        predictionMode_block = invzigzag(predictionMode1d, 1, 1);

                        % Apply differential decoding
                        [predictionMode_block, previous_prediction_mode_block] = ...
                            diffDecoding_block(predictionMode_block, 'modes', previous_prediction_mode_block);

                        % Assign the mode to all positions in the 2x2 block
                        decodedPredicitonModes2d(row:row+1, col:col+1) = predictionMode_block;

                        % Set vbs_matrix to 0 (large block)
                        reconstructed_vbs_matrix(row:row+1, col:col+1) = 0;

                    else
                        error('Invalid prefix encountered during decoding!');
                    end
                end
            end
        else
            predictedModeRevEGC = exp_golomb_decode(encodedPredicitonModes1d); 
            decodedPredicitonModes2d = invzigzag(predictedModeRevEGC, predwidth, predheight);
            decodedPredicitonModes2d = diffDecoding(decodedPredicitonModes2d,'modes');
        end
    end
    
    reconstructedResiduals2d = zeros(reswidth, resheight);
    decodedResidues2d = zeros(reswidth, resheight);
    idx = 1;  % Index for traversing encodedResidues1d
    temp = reswidth;
    reswidth = resheight;
    resheight = temp;

    % Assuming dct_blockSize is known (should be consistent with quantization_entropy function)
    if VBSEnable
        vbs_matrix = reconstructed_vbs_matrix;
        [vbsRows, vbsCols] = size(vbs_matrix);

        for blockY = 1:2:vbsRows
            for blockX = 1:2:vbsCols
                rowOffset = (blockY -1) * dct_blockSize + 1;
                colOffset = (blockX -1) * dct_blockSize + 1;
                actualBlockHeight = min(2 * dct_blockSize, resheight - rowOffset +1);
                actualBlockWidth = min(2 * dct_blockSize, reswidth - colOffset +1);
                vbs_block = vbs_matrix(blockY:blockY+1, blockX:blockX+1);

                if all(vbs_block(:) == 0)
                    % Process as a large block
                    % Skip any -1 separators
                    while idx <= length(encodedResidues1d) && encodedResidues1d(idx) == -1
                        idx = idx + 1;
                    end
                    % Collect encodedResidues1d for this block until next -1 or end
                    encodedResidues_block_wth_QP = [];
                    while idx <= length(encodedResidues1d) && encodedResidues1d(idx) ~= -1
                        encodedResidues_block_wth_QP = [encodedResidues_block_wth_QP, encodedResidues1d(idx)];
                        idx = idx +1;
                    end
                    baseQP = encodedResidues_block_wth_QP(1);
                    encodedResidues_block = encodedResidues_block_wth_QP(2:end);
                    % Decode this block
                    residuesRevEGC = exp_golomb_decode(encodedResidues_block);
                    numElements = actualBlockHeight * actualBlockWidth;
                    decodedResidues1d = rle_decode(residuesRevEGC, numElements);
                    if isempty(decodedResidues1d)
                        break;
                    end
                    decodedResidues2d_block = invzigzag(decodedResidues1d, actualBlockHeight, actualBlockWidth);

                     % Create the quantization matrix
                    Q = createQMatrix(size(decodedResidues2d_block), baseQP);
                    
                    % Inverse quantization (element-wise multiplication)
                    dequantizedSubBlock = decodedResidues2d_block .* Q;
                    
                    % Apply inverse DCT to the dequantized sub-block
                    idctSubBlock = idct2(dequantizedSubBlock);

                    % Place into decodedResidues2d
                    decodedResidues2d(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = decodedResidues2d_block;
                     reconstructedResiduals2d(rowOffset:rowOffset+actualBlockHeight-1, colOffset:colOffset+actualBlockWidth-1) = idctSubBlock;
                else
                    % Process each sub-block separately
                    for subBlockY = 0:1
                        for subBlockX = 0:1
                            subRowOffset = rowOffset + subBlockY * dct_blockSize;
                            subColOffset = colOffset + subBlockX * dct_blockSize;
                            if subRowOffset > resheight || subColOffset > reswidth
                                continue;
                            end
                            actualSubBlockHeight = min(dct_blockSize, resheight - subRowOffset +1);
                            actualSubBlockWidth = min(dct_blockSize, reswidth - subColOffset +1);
                            % Skip any -1 separators
                            while idx <= length(encodedResidues1d) && encodedResidues1d(idx) == -1
                                idx = idx +1;
                            end
                            % Collect encodedResidues1d for this block until next -1 or end
                            encodedResidues_block_wth_QP = [];
                            while idx <= length(encodedResidues1d) && encodedResidues1d(idx) ~= -1
                                encodedResidues_block_wth_QP = [encodedResidues_block_wth_QP, encodedResidues1d(idx)];
                                idx = idx +1;
                            end
                            if isempty(encodedResidues_block_wth_QP)
                                break;
                            end
                            baseQP = encodedResidues_block_wth_QP(1);
                            encodedResidues_block = encodedResidues_block_wth_QP(2:end);
                            % Decode
                            residuesRevEGC = exp_golomb_decode(encodedResidues_block);
                            numElements = actualSubBlockHeight * actualSubBlockWidth;
                            decodedResidues1d = rle_decode(residuesRevEGC, numElements);
                            if isempty(decodedResidues1d)
                                break;
                            end
                            decodedResidues2d_block = invzigzag(decodedResidues1d, actualSubBlockHeight, actualSubBlockWidth);
                            % Place into decodedResidues2d
                            % Create the quantization matrix
                            Q = createQMatrix(size(decodedResidues2d_block), baseQP);
                            
                            % Inverse quantization (element-wise multiplication)
                            dequantizedSubBlock = decodedResidues2d_block .* Q;
                            
                            % Apply inverse DCT to the dequantized sub-block
                            idctSubBlock = idct2(dequantizedSubBlock);
                            decodedResidues2d(subRowOffset:subRowOffset+actualSubBlockHeight-1, subColOffset:subColOffset+actualSubBlockWidth-1) = decodedResidues2d_block;
                             reconstructedResiduals2d(subRowOffset:subRowOffset+actualSubBlockHeight-1, subColOffset:subColOffset+actualSubBlockWidth-1) = idctSubBlock;
                        end
                    end
                end
            end
        end
    else
        % Decoding residues without VBS
        for row = 1:dct_blockSize:resheight
            for col = 1:dct_blockSize:reswidth
                actualBlockHeight = min(dct_blockSize, resheight - row + 1);
                actualBlockWidth = min(dct_blockSize, reswidth - col + 1);
                % Skip any -1 separators
                while idx <= length(encodedResidues1d) && encodedResidues1d(idx) == -1
                    idx = idx + 1;
                end

                encodedResidues_block_wth_QP = [];
                while idx <= length(encodedResidues1d) && encodedResidues1d(idx) ~= -1
                    encodedResidues_block_wth_QP = [encodedResidues_block_wth_QP, encodedResidues1d(idx)];
                    idx = idx +1;
                end
                baseQP = encodedResidues_block_wth_QP(1);
                encodedResidues_block = encodedResidues_block_wth_QP(2:end);


                % Collect encodedResidues1d for this block until next -1 or end
                % encodedResidues_block = [];
                % while idx <= length(encodedResidues1d) && encodedResidues1d(idx) ~= -1
                %     encodedResidues_block = [encodedResidues_block, encodedResidues1d(idx)];
                %     idx = idx +1;
                % end
                % Decode
                residuesRevEGC = exp_golomb_decode(encodedResidues_block);
                numElements = actualBlockHeight * actualBlockWidth;
                decodedResidues1d = rle_decode(residuesRevEGC, numElements);
                if isempty(decodedResidues1d)
                    break;
                end
                decodedResidues2d_block = invzigzag(decodedResidues1d, actualBlockHeight, actualBlockWidth);
                 % Create the quantization matrix
                Q = createQMatrix(size(decodedResidues2d_block), baseQP);
                
                % Inverse quantization (element-wise multiplication)
                dequantizedSubBlock = decodedResidues2d_block .* Q;
                
                % Apply inverse DCT to the dequantized sub-block
                idctSubBlock = idct2(dequantizedSubBlock);
                % Place into decodedResidues2d
                decodedResidues2d(row:row+actualBlockHeight-1, col:col+actualBlockWidth-1) = decodedResidues2d_block;
                reconstructedResiduals2d(row:row+actualBlockHeight-1, col:col+actualBlockWidth-1) = idctSubBlock;
            end
        end
    end

end


function matrix_3d = reshape_2d_to_3d(reshaped_2d, rows, cols, depth)
    % Reshape the 2D matrix back to the original 3D matrix
    % Each row block in the original 3D matrix is represented in flattened form in reshaped_2d
    matrix_3d = reshape(reshaped_2d, cols, rows, depth);
end



function matrix = invzigzag(zigzag_order, rows, cols)
    % Inverse Zigzag function to reorder a 1D array back into a 2D matrix
    matrix = zeros(rows, cols);
    index = 1;
    for s = 1:(rows + cols - 1)
        % Even diagonals: top-right to bottom-left
        for i = max(1, s - cols + 1):min(rows, s)
            j = s - i + 1;
            matrix(i, j) = zigzag_order(index);
            index = index + 1;
        end
    end
end



% function encoded = exp_golomb_encode(data)
%     % Encode data using Exponential-Golomb coding based on specific rules
%     % and return a 1D cell array of binary strings
%     encoded = cell(1, length(data));  % Initialize a cell array to store binary strings
%     for i = 1:length(data)
%         x = data(i);
%         if x > 0
%             value = 2 * x;  % Positive x to odd integer (2x-1)
%         else
%             value = (-2 * x)+1;     % Non-positive x to even integer (-2x)
%         end
%         bin_code = dec2bin(value);  % Convert to binary
%         leading_zeros = floor(log2(value)) + 1;  % Compute the number of leading zeros
%         code = [repmat('0', 1, leading_zeros - 1), '1', bin_code(2:end)];  % Assemble the code
%         encoded{i} = code;  % Store binary code in cell array
%     end
% end

function data = exp_golomb_decode(encoded)
    % Decode a 1D array of binary values (0s and 1s) encoded using Exponential-Golomb coding
    data = [];  % Initialize an empty array to store the decoded values
    idx = 1;
    
    while idx <= length(encoded)
        % Count the number of leading zeros
        leading_zeros = 0;
        while idx <= length(encoded) && encoded(idx) == 0
            leading_zeros = leading_zeros + 1;
            idx = idx + 1;
        end
        
        % Read the remainder of the binary value
        if idx <= length(encoded) && encoded(idx) == 1
            idx = idx + 1;
            bin_code_length = leading_zeros;
            bin_code = [1, encoded(idx : idx + bin_code_length - 1)];
            idx = idx + bin_code_length;
            
            % Convert binary code to decimal value
            value = bin2dec(char(bin_code + '0'));
            
            % Map value back to original data
            if mod(value, 2) == 0
                x = (value / 2);
            else
                x = - (value - 1) / 2;
            end
            
            data = [data, x];  % Append the decoded value to the output array
        end
    end
end


function decoded = rle_decode(encoded,datalength)
    % Decode data using a custom Run-Length Encoding (RLE) scheme
    % Return a 1D array containing the decoded values of shape 1x(blockSize * blockSize)
    decoded = [];
    i = 1;
    
    while i <= length(encoded)
        if encoded(i) < 0
            % Negative value indicates a run of non-zero values
            non_zero_count = -encoded(i);
            i = i + 1;
            decoded = [decoded, encoded(i:i + non_zero_count - 1)];
            i = i + non_zero_count;
        elseif encoded(i) > 0
            % Positive value indicates a run of zeros
            zero_count = encoded(i);
            
            decoded = [decoded, zeros(1, zero_count)];
            i = i + 1;
        else
            % Zero value in encoded data indicates an actual zero in the original data
            decoded = [decoded, zeros(1, datalength-length(decoded))];
            i = length(encoded) + 1;
        end
    end
end


function encoded = rle_encode(data)
    % Encode data using a custom Run-Length Encoding (RLE) scheme
    % Return a 1D array containing the encoded values
    encoded = [];
    
    i = 1;
    while i <= length(data)
        if data(i) ~= 0
            % Count the number of consecutive non-zero values
            non_zero_count = 0;
            start_idx = i;
            while i <= length(data) && data(i) ~= 0
                non_zero_count = non_zero_count + 1;
                i = i + 1;
            end
            % Encode the non-zero run: prefix with -count followed by the values
            encoded = [encoded, -non_zero_count, data(start_idx:start_idx + non_zero_count - 1)];
        else
            % Count the number of consecutive zero values
            zero_count = 0;
            while i <= length(data) && data(i) == 0
                zero_count = zero_count + 1;
                i = i + 1;
            end

            if i == length(data)+1
                zero_count = 0;
            end
            % Encode the zero run: prefix with +count
            encoded = [encoded, zero_count];
        end
    end
    
    % Add a trailing 0 if the rest of the elements are zeros
end

function Q = createQMatrix(blockSize, QP)
    if numel(blockSize) > 1
        rows = blockSize(1);
        cols = blockSize(2);
    else
        rows = blockSize;
        cols = blockSize;
    end
    
    Q = zeros(rows, cols);
    for x = 1:rows
        for y = 1:cols
            if (x + y <= rows)
                Q(x,y) = 2^QP;
            elseif (x + y == rows + 1)
                Q(x,y) = 2^(QP+1);
            else
                Q(x,y) = 2^(QP+2);
            end
        end
    end
end

% testdata = [-31, 9, -4, 8, 1, -3, 4, 4, 2, 4, 0, 4, 0, 0, -4, 0, 0, 1, 0, 0]
% encoded = rle_encode(testdata)
% decoded = rle_decode(encoded,20)